#ifndef MC_PARTICLE_BUFFER_INCLUDE
#define MC_PARTICLE_BUFFER_INCLUDE

#include "MC_Processor_Info.hh"
#include "utilsMpi.hh"
#include <map>

// forward declarations
class MonteCarlo;

class mcp_test_done_class
{
 public:
    mcp_test_done_class() { Zero_Out(); }
    int local_sent;
    int local_recv;
    int BlockingSum;

    void Get_Local_Gains_And_Losses(MonteCarlo *mcco, int64_t sent_recv[2]);
    void Post_Recv();
    void Zero_Out();
    void Reduce_Num_Sent_And_Recv(int64_t buf_sum[2]);
    void Free_Memory();
    bool ThisProcessorCommunicates(int rank = -1);
};


class MC_Particle_Buffer
{
 private:

    MonteCarlo *mcco;
    mcp_test_done_class          test_done;

    void Instantiate();
    void Initialize_Map();
    bool Trivially_Done();

 public:
    // non-master threads place full buffers here for master thread to send
    // std::list<particle_buffer_base_type> thread_send_buffer_queue;

    std::map<int, int>    processor_buffer_map; // Map processors to buffers. buffer_index = processor_buffer_map[processor]
    int  num_buffers;         // Number of particle buffers
    int  buffer_size;         // Buffer size to be sent.

    MC_Particle_Buffer(MonteCarlo *mcco_, size_t bufferSize_);       // constructor
    void Initialize();
    int  Choose_Buffer(int neighbor_rank);
    bool Test_Done_New();
    bool Allreduce_ParticleCounts();

    void Free_Memory();
private:
    MC_Particle_Buffer( const MC_Particle_Buffer& );                    // disable copy constructor
    MC_Particle_Buffer& operator=( const MC_Particle_Buffer& tmp );     // disable assignment operator
};

/*
  mcco->particle_buffer->test_done.local_sent
                                  .local_recv
                                  .BlockingSum;



                         num_buffers
                         task[]
                               .send_buffer[].processor
                                             .num_particles
                                             .int_index, float_index, char_index, length
                                             .int_data, float_data, char_data
                                             .request_list
                               .recv_buffer[].processor
                                             .num_particles
                                             .int_index, float_index, char_index, length
                                             .int_data, float_data, char_data
                                             .request_list
                         processor_buffer_map[]

 */
#endif
