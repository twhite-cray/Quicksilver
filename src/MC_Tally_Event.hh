#pragma once
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

