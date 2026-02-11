--  ============================================================================
--  GEM.Random_Descent - Stochastic Optimization via Randomized Parameter Search
--  ============================================================================
--
--  PURPOSE:
--    Implements a randomized descent optimization algorithm for parameter
--    search in high-dimensional space. Uses probabilistic perturbations to
--    explore parameter space and escape local minima.
--
--  ALGORITHM:
--    Markov Chain Monte Carlo (MCMC) style random walk:
--    1. Select a random parameter from the set
--    2. Apply logarithmic or linear perturbation (spread-controlled)
--    3. Accept or reject based on fitness improvement (done by caller)
--    4. Optionally flip sign to explore negative parameter space
--
--  KEY FEATURES:
--    - Logarithmic perturbations: Better exploration across orders of magnitude
--    - Relative vs. Absolute adjustment modes (REL environment variable)
--    - Sign flipping: FLIP parameter allows exploring negative values
--    - Harmonic randomization: Mutation of frequency harmonics
--    - Fixed parameter support: Some parameters locked during optimization
--
--  RANDOM NUMBER GENERATORS:
--    - G (Float_Random): Continuous values for parameter adjustments
--    - D (Discrete_Random): Selects which parameter to perturb
--    - H (Discrete_Random): Selects harmonic frequencies
--
--  NOTE ON REPRODUCIBILITY:
--    Currently uses default random seed (time-based). For reproducible
--    testing, generators could be seeded with Reset(Gen, Initiator) using
--    a fixed value from SEED environment variable. This would allow
--    deterministic floating-point results across runs.
--
--  CONFIGURATION:
--    - FLIP: Probability of sign flip (0.0 = never, 1.0 = always)
--    - FIX: Lock harmonics during optimization
--    - LOGH: Use logarithmic distribution for harmonic selection
--    - REL: Relative (multiplicative) vs absolute (additive) adjustments
--    - RESET: Reset RNG to time-based seed on startup
--
--  ============================================================================

with Ada.Numerics.Float_Random;
with Ada.Numerics.Long_Elementary_Functions;
with Text_IO;
with Ada.Long_Float_Text_IO;
with Ada.Numerics.Discrete_Random;

package body GEM.Random_Descent is
   package FR renames Ada.Numerics.Float_Random;
   package LEF renames Ada.Numerics.Long_Elementary_Functions;

   subtype Set_Index is Positive range 1 .. Set_Range;
   package DR is new Ada.Numerics.Discrete_Random (Set_Index);
   Flip_Value : constant Long_Float := GEM.Getenv ("FLIP", 0.0);
   Fix_Harm : constant Boolean := GEM.Getenv ("FIX", False);
   Log_Harm : constant Boolean := GEM.Getenv ("LOGH", False);
   Relative : constant Boolean := GEM.Getenv ("REL", True);

   subtype Harmonic_Index is Positive range 2 .. Harmonic_Range;
   package HR is new Ada.Numerics.Discrete_Random (Harmonic_Index);

   --  Random number generators (time-seeded by default)
   --  TODO: For reproducible testing, could add seeding via SEED env var:
   --    if GEM.Getenv("SEED", 0) /= 0 then
   --       DR.Reset(D, DR.Initiator(GEM.Getenv("SEED", 0)));
   --       FR.Reset(G, FR.Initiator(GEM.Getenv("SEED", 0)));
   --       HR.Reset(H, HR.Initiator(GEM.Getenv("SEED", 0)));
   --    end if;
   --  This would enable deterministic floating-point results for regression testing.
   D : DR.Generator; -- Discrete for selecting from a set of params
   G : FR.Generator; -- Floating point for values
   H : HR.Generator; -- Floating point for values

   --  =========================================================================
   --  Markov: Apply random perturbation to parameter array
   --
   --  STRATEGY:
   --    Randomly select one parameter and perturb it using logarithmic or
   --    linear scaling. Skips zero values and fixed parameters. Uses
   --    recursion to retry if selected parameter cannot be modified.
   --
   --  PERTURBATION MODES:
   --    - Relative (REL=True): Multiplicative adjustment scaled by log(random)
   --      New_Value = Old_Value * (1 + Spread * log(Ran) * Sign)
   --    - Absolute (REL=False): Additive adjustment from calibration baseline
   --      New_Value = Cal_Value + Spread * Cal_Value * Ran * Sign
   --
   --  PARAMETERS:
   --    Set: Parameter array to modify (in/out)
   --    Ref: Backup copy before modification (out)
   --    Spread: Controls perturbation magnitude (typically 0.0001 to 0.1)
   --    Cal: Calibration baseline values for absolute mode
   --  =========================================================================
   procedure Markov
     (Set : in out LF_Array; Ref : out LF_Array; Spread : in Long_Float;
      Cal : in LF_Array)
   is
      Ran : Long_Float := Long_Float (FR.Random (G));
      Sign : Long_Float :=
        Long_Float (FR.Random (G) - 0.5); -- Is step + or - ?
      Adjust : Long_Float;
      I : Set_Index := DR.Random (D);
   begin
      Ref := Set; -- keep in case of error?
      --  TODO: Can remove - Original gradient descent idea was to update all
      --  parameters simultaneously (full gradient). Random descent instead
      --  updates one parameter at a time (line 51) for simpler Markov chain.
      --  Full gradient would require computing sensitivity of all parameters.
      --for I in Set'Range loop -- a real gradient descent would do all at once
      if Set (I) = 0.0 then
         if Flip_Value < 0.0 then
            Set (I) := Spread * Long_Float'Copy_Sign (LEF.Log (Ran), Sign);
         else
            Markov (Set, Ref, Spread, Cal);  -- recurse, pick another
         end if;
      elsif Fixed (Set (I)) then
         Markov (Set, Ref, Spread, Cal);  -- recurse, pick another
      else
         if Relative then
            Adjust :=
              1.0 + Spread * Long_Float'Copy_Sign (LEF.Log (Ran), Sign);
            Set (I) := Set (I) * Adjust;
            if Ran < Flip_Value then
               Set (I) := -Set (I);
            end if;
         else
            Adjust := Spread * Cal (I) * Long_Float'Copy_Sign (Ran, Sign);
            Set (I) := Cal (I) + Adjust;
         end if;
      end if;
      --  TODO: Can remove - End of gradient descent loop (see above)
      --end loop;
   end Markov;

   --  Overloaded version: Apply random perturbation to single scalar value
   --  Same algorithm as array version but operates on one value directly
   procedure Markov
     (Value : in out Long_Float; Ref : out Long_Float; Spread : in Long_Float;
      Cal : in Long_Float)
   is
      Ran : Long_Float := Long_Float (FR.Random (G));
      Sign : Long_Float :=
        Long_Float (FR.Random (G) - 0.5); -- Is step + or - ?
      Adjust : Long_Float;
   begin
      Ref := Value;
      if Relative then
         Adjust := 1.0 + Spread * Long_Float'Copy_Sign (LEF.Log (Ran), Sign);
         Value := Value * Adjust;
         if Ran < Flip_Value then
            Value := -Value;
         end if;
      else
         Adjust := Spread * Cal * Long_Float'Copy_Sign (Ran, Sign);
         Value := Cal + Adjust;
      end if;
   end Markov;

   --  Debug utility: Print all parameters in array
   procedure Dump (Set : in LF_Array) is
   begin
      for I in Set'Range loop
         Ada.Long_Float_Text_IO.Put (Set (I), Fore => 4, Aft => 11, Exp => 0);
         Text_IO.Put_Line (", ");
      end loop;
   end Dump;

   --  =========================================================================
   --  Random_Harmonic: Randomize frequency harmonic indices
   --
   --  Used to explore different harmonic combinations in Fourier-like
   --  decomposition. Can use uniform (HR.Random) or logarithmic (Log)
   --  distribution to favor lower harmonics.
   --
   --  FIX_HARM: Lock harmonics (no randomization)
   --  LOG_HARM: Use logarithmic weighting (favors lower harmonics)
   --  =========================================================================
   procedure Random_Harmonic (Index : in out Positive; Ref : out Positive) is
      Ran : Long_Float := Long_Float (FR.Random (G));
      FV : Long_Float := abs Flip_Value;
   begin
      Ref := Index;
      if Fix_Harm then
         null;
      elsif Ran < FV then
         if Log_Harm then
            Index :=
              2 -
              Integer
                (Long_Float (Harmonic_Range) *
                 LEF.Log (Long_Float (FR.Random (G))));
         else
            Index := HR.Random (H);
         end if;
      end if;
   end Random_Harmonic;

   --  Overloaded version: Randomize array of harmonics with duplicate checking
   --  Ensures no duplicate harmonic indices (would cause singular matrix in
   --  multivariate regression). Retries if duplicate detected.
   procedure Random_Harmonic (Index : in out Ns; Ref : out Ns) is
      Test : Positive;
      Reference : Ns := Index;
   begin
      for I in Index'Range loop
         Random_Harmonic (Index (I), Test);
         for J in Reference'Range loop
         -- Make sure no duplicates in Harmonics otherwise singularity in sol'n
            if Index (I) = Reference (J) then
               Index (I) := Test;  -- keep the old value
               exit;
            end if;
         end loop;
      end loop;
      Ref := Reference;
   end Random_Harmonic;

   --  Reset all random number generators to time-based seed
   --  Used with RESET environment variable for fresh start
   procedure Reset is
   begin
      DR.Reset (D);
      FR.Reset (G);
      HR.Reset (H);
   end Reset;

   --  Generate small random noise for mixing into forcing signal
   --  Adds jitter scaled by 1/1000 to avoid numerical stagnation
   function Small_Random (Last : in Long_Float) return Long_Float is
   begin
      return Last + (Long_Float (FR.Random (G)) - 0.5) / 1_000.0;
   end Small_Random;

begin
   --  Optional reset on package initialization
   if GEM.Getenv ("RESET", False) then
      Reset;
   end if;
end GEM.Random_Descent;
