#ifndef TALLIES_HH
#define TALLIES_HH

#include "portability.hh"
#include "QS_Vector.hh"
#include <stdlib.h>
#include <string>
#include <vector>
#include <cinttypes>
#include "NuclearData.hh"
#include "MonteCarlo.hh"
#include "utils.hh"
#include "macros.hh"
#include "BulkStorage.hh"
#include "DeclareMacro.hh"
#include "EnergySpectrum.hh"

typedef unsigned long long int uint64_cu;

class Fluence;

struct MC_Tally_Event
{
   enum Enum
   {
      Collision,
      Facet_Crossing_Transit_Exit,
      Census,
      Facet_Crossing_Tracking_Error,
      Facet_Crossing_Escape,
      Facet_Crossing_Reflection,
      Facet_Crossing_Communication
   };
};

class Balance
{
   public:

   uint64_cu _absorb;      // Number of particles absorbed
   uint64_cu _census;      // Number of particles that enter census
   uint64_cu _escape;      // Number of particles that escape
   uint64_cu _collision;   // Number of collosions
   uint64_cu _end;         // Number of particles at end of cycle
   uint64_cu _fission;     // Number of fission events
   uint64_cu _produce;     // Number of particles created by collisions
   uint64_cu _scatter;     // Number of scatters
   uint64_cu _start;       // Number of particles at beginning of cycle
   uint64_cu _source;      // Number of particles sourced in
   uint64_cu _rr;          // Number of particles Russian Rouletted in population control
   uint64_cu _split;       // Number of particles split in population control
   uint64_cu _numSegments; // Number of segements

   Balance() :
      _absorb(0), _census(0), _escape(0), _collision(0), _end(0), _fission(0), _produce(0), _scatter(0), _start(0),
      _source(0), _rr(0), _split(0), _numSegments(0) { }

   ~Balance() {}

    void PrintHeader()
    {
        Print0("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s",
               "start", "source", "rr", "split", 
               "absorb", "scatter", "fission", "produce", 
               "collisn", "escape", "census", 
               "num_seg");
    }

    void Print()
    {
        Print0("%12" PRIu64 " %12" PRIu64 " %12" PRIu64 " %12" PRIu64 " %12" PRIu64 " %12" PRIu64 " %12" PRIu64 " %12" PRIu64 " %12" PRIu64 " %12" PRIu64 " %12" PRIu64 " %12" PRIu64 "",
               _start, _source,_rr, _split,
               _absorb, _scatter, _fission, _produce, 
               _collision, _escape, _census, 
               _numSegments);
    }

   void Reset()
   {
      _absorb = _census = _escape = _collision = _end = _fission = _produce = _scatter = _start = _source =
          _rr = _split = _numSegments = 0;
   }

   void Add(Balance &bal)
   {
      _absorb      += bal._absorb;
      _census      += bal._census;
      _escape      += bal._escape;
      _collision   += bal._collision;
      _end         += bal._end;
      _fission     += bal._fission;
      _produce     += bal._produce;
      _scatter     += bal._scatter;
      _start       += bal._start;
      _source      += bal._source;
      _rr          += bal._rr;
      _split       += bal._split;
      _numSegments += bal._numSegments;
   }
};

class ScalarFluxCell
{
   public:
   double* _group;
   int _size;

   ScalarFluxCell() : _group(0), _size(0) {}

   ScalarFluxCell(double* storage, int size)
   :  _group(storage),
      _size(size)
   {
      for( int i = 0; i < _size; i++ )
      {
         _group[i] = 0.0;
      }
   }
   
   ~ScalarFluxCell() {}

   int size() const {return _size;}
};


class CellTallyTask
{
  public:
    qs_vector<double> _cell;

    CellTallyTask() : _cell() {}

    CellTallyTask(MC_Domain* domain)
    {
        if( _cell.capacity() == 0 )
        {
            _cell.reserve( domain->cell_state.size(), VAR_MEM );
        }

        _cell.Open();
        for (int cellIndex = 0; cellIndex < domain->cell_state.size(); cellIndex++)
        {
            _cell.push_back(0.0);
        }
        _cell.Close();
    }

    void Add(CellTallyTask &cellTallyTask)
    {
        for (int cellIndex = 0; cellIndex < _cell.size(); cellIndex++)
        {
            _cell[cellIndex] += cellTallyTask._cell[cellIndex];
        }
    }

    void Reset() 
    {
        for (int cellIndex = 0; cellIndex < _cell.size(); cellIndex++)
        {
            _cell[cellIndex] = 0.0;
        }     
    }

    ~CellTallyTask() {}
};

class ScalarFluxTask
{
   public:
   qs_vector<ScalarFluxCell> _cell;
   BulkStorage<double> _scalarFluxCellStorage;

   ScalarFluxTask() : _cell() {}

   ScalarFluxTask(MC_Domain* domain, int numGroups)
   {
       if( _cell.capacity() == 0 )
       {
          _cell.reserve( domain->cell_state.size(), VAR_MEM );
          _scalarFluxCellStorage.setCapacity(domain->cell_state.size()*numGroups, VAR_MEM);
       }

       _cell.Open();
      for (int cellIndex = 0; cellIndex < domain->cell_state.size(); cellIndex++)
      {
         double* tmp = _scalarFluxCellStorage.getBlock(numGroups);
         _cell.push_back(ScalarFluxCell(tmp, numGroups));
      }
      _cell.Close();
   }

   void Add(ScalarFluxTask &scalarFluxTask)
   {
      unsigned int numGroups = _cell[0].size();
      for (int cellIndex = 0; cellIndex < _cell.size(); cellIndex++)
      {
         for (int groupIndex = 0; groupIndex < numGroups; groupIndex++)
         {
             _cell[cellIndex]._group[groupIndex] += scalarFluxTask._cell[cellIndex]._group[groupIndex];
         }
      }
   }

   void Reset() 
   {
      unsigned int numGroups = _cell[0].size();
      for (int cellIndex = 0; cellIndex < _cell.size(); cellIndex++)
      {
         for (int groupIndex = 0; groupIndex < numGroups; groupIndex++)
         {
            _cell[cellIndex]._group[groupIndex] = 0.0;
         }
      }     
   }

   ~ScalarFluxTask() {}
};

class CellTallyDomain
{
   public:
   qs_vector<CellTallyTask> _task;

   CellTallyDomain() : _task() {}

   CellTallyDomain(MC_Domain* domain)
   {
      // Assume OMP_NUM_THREADS tasks
      if( _task.capacity() == 0 )
      {
        _task.reserve(1, VAR_MEM);
      }
      _task.Open();
      _task.push_back(CellTallyTask(domain));
      _task.Close();
   }

   ~CellTallyDomain() {}
};

class ScalarFluxDomain
{
   public:
   qs_vector<ScalarFluxTask> _task;

   ScalarFluxDomain() : _task() {}

   ScalarFluxDomain(MC_Domain* domain, int numGroups)
   {
      // Assume OMP_NUM_THREADS tasks
      if( _task.capacity() == 0 )
      {
        _task.reserve(1, VAR_MEM);
      }
      _task.Open();
      _task.push_back(ScalarFluxTask(domain, numGroups));
      _task.Close();
   }

   ~ScalarFluxDomain() {}
};

class FluenceDomain
{
  public:
    FluenceDomain( int numCells) : _cell(numCells, 0.0)
    {}

    void addCell( int index, double value ){ _cell[index] += value;}
    double getCell( int index ){ return _cell[index]; }
    int size(){ return _cell.size(); }
    
  private:
    std::vector<double> _cell;
};

class Fluence
{
  public:
    Fluence() {};
    ~Fluence()
    {
        for( int i = 0; i < _domain.size(); i++ )
        {
            if( _domain[i] != NULL )
                delete _domain[i];
        }
    }

    void compute(int domain, ScalarFluxDomain &scalarFluxDomain);

    std::vector<FluenceDomain*> _domain;
};

class Tallies
{
  public:
    Balance                     _balanceCumulative;
    qs_vector<Balance>          _balanceTask;
    qs_vector<ScalarFluxDomain> _scalarFluxDomain;
    qs_vector<CellTallyDomain>  _cellTallyDomain;
    Fluence                     _fluence;
    EnergySpectrum              _spectrum;
    
    Tallies( std::string spectrumName, int spectrumSize )
      : _balanceCumulative(), _balanceTask(),
        _scalarFluxDomain(), _spectrum(spectrumName, spectrumSize)
    {
    }

    ~Tallies() {}

    void InitializeTallies( MonteCarlo *monteCarlo ); 

    void CycleInitialize(MonteCarlo* monteCarlo);

    void CycleFinalize(MonteCarlo *mcco);
    void PrintSummary(MonteCarlo *mcco);

    void TallyCellValue(double value, int domain, int cell)
    {
        ATOMIC_ADD( _cellTallyDomain[domain]._task[0]._cell[cell], value );
    }

    double ScalarFluxSum(MonteCarlo *mcco);
};

#endif
