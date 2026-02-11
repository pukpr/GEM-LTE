--  ============================================================================
--  GEM.LTE - Laplace's Tidal Equation Package Body
--  ============================================================================
--
--  PURPOSE:
--    Implements Laplace's Tidal Equation (LTE) framework for modeling
--    ocean-atmosphere response to long-period gravitational forcing.
--    Provides Doodson number calculations and initialization of tidal
--    constituent sets used throughout GEM-LTE climate analysis.
--
--  SCIENTIFIC BACKGROUND:
--    Doodson numbers encode the astronomical arguments (s, h, p, N, etc.)
--    that determine tidal periods. Each constituent represents a specific
--    combination of lunar/solar forcing frequencies. The framework supports:
--    - Standard long-period tides (18.6yr, 8.85yr, etc.)
--    - Annual/semi-annual harmonics
--    - Quasi-Biennial Oscillation (QBO) periods
--
--  KEY FUNCTIONS:
--    - Doodson: Converts Doodson arguments to frequency (cycles/year)
--    - Year_Adjustment: Dynamically updates all periods when year length changes
--
--  INITIALIZATION:
--    Package elaboration computes all tidal periods and initializes:
--    - LP/LPAP: Main long-period constituents
--    - LP_Set/LPAP_Set: Extended constituent set
--    - LP_RSet/LPAP_RSet: Reduced constituent set
--    - LP_QSet/LPAP_QSet: QBO-focused set
--    - LP_Annual/LPAP_Annual: Annual harmonics
--    - QBO/QBOAP: QBO-specific periods
--
--  ============================================================================

package body GEM.LTE is

   Year_in_Days : constant := 365.242_248_4;  -- 365.241237718675000;
   Year_Correction : Long_Float := GEM.Getenv ("YEAR", 0.0);
   --  TODO: Can remove - dLOD_Mod was experimental scaling factor for day-length
   --  variations. Standard Doodson calculation proved sufficient.
   --dLOD_Mod : Long_Float := GEM.Getenv("dLOD_Mod", 1.0);
   Year_Dynamic_Correction : Long_Float := 0.0;

   function Year_Length return Long_Float is
   begin
      return Year_in_Days + Year_Correction + Year_Dynamic_Correction;
   end Year_Length;

   --  Compute tidal frequency from Doodson arguments
   --  Formula combines tropical year, lunar nodal, perigee, and anomalistic periods
   function Doodson (I : in Integer; D : in Doodson_List) return Long_Float is
   begin
      --  TODO: Can remove - Conditional dLOD scaling was experimental approach
      --  to handle day-length variations. Standard formula works for all cases.
      --if D(I).h = 1 and dLOD_Mod /= 1.0 then
      --   return dLOD_Mod * Year_Length;
      --else
      return
        1.0 /
        (Long_Float (D (I).s) / Tropical +
         (Long_Float (D (I).h) * D (I).Period) / Year_Length +
         Long_Float (D (I).p) / p + Long_Float (D (I).N) / N);
      --end if;
   end Doodson;

   --  Dynamically adjust all tidal periods when year length changes
   --  Used during optimization to refine astronomical constants
   procedure Year_Adjustment (Value : in Long_Float; List : in out Periods) is
   begin
      if Value /= 0.0 then
         Year_Dynamic_Correction := Value;
         for I in List'Range loop
            List (I) := Doodson (I, Doodson_Args);
         end loop;
      end if;
   end Year_Adjustment;

begin

   --  Package elaboration: Initialize all tidal constituent period arrays
   --  Converts Doodson arguments to frequencies for use in LTE solver

   for I in Doodson_Args'Range loop
      Doodson_Args (I).Period := Doodson (I, Doodson_Args);
      LP (I) := Doodson_Args (I).Period;
      LPAP (I) := (0.01, 1.0);
   end loop;

   for I in Doodson_Set'Range loop
      Doodson_Set (I).Period := Doodson (I, Doodson_Set);
      LP_Set (I) := Doodson_Set (I).Period;
      LPAP_Set (I) := (0.01, 1.0);
   end loop;

   for I in Doodson_RSet'Range loop
      Doodson_RSet (I).Period := Doodson (I, Doodson_RSet);
      LP_RSet (I) := Doodson_RSet (I).Period;
      LPAP_RSet (I) := (0.01, 1.0);
   end loop;

   for I in Doodson_QSet'Range loop
      Doodson_QSet (I).Period := Doodson (I, Doodson_QSet);
      LP_QSet (I) := Doodson_QSet (I).Period;
      LPAP_QSet (I) := (0.01, 1.0);
   end loop;

   for I in Annual_Set'Range loop
      Annual_Set (I).Period := Doodson (I, Annual_Set);
      LP_Annual (I) := 1.0 / Annual_Set (I).Period;
      LPAP_Annual (I) := (0.01, 1.0);
   end loop;

   for I in QBO_Args'Range loop
      QBO_Args (I).Period := Doodson (I, QBO_Args);
      QBO (I) := QBO_Args (I).Period;
      QBOAP (I) := (0.01, 1.0);
   end loop;

end GEM.LTE;
