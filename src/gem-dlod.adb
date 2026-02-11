--  ============================================================================
--  GEM.dLOD - Day Length of Day (dLOD) Tidal Analysis
--  ============================================================================
--
--  PURPOSE:
--    Analyzes variations in Earth's rotation rate (day length) caused by
--    gravitational tidal forces. Performs multivariate regression to extract
--    amplitude and phase of long-period tidal constituents from observed dLOD
--    data. This provides validation that GEM-LTE can correctly model known
--    astronomical forcing.
--
--  SCIENTIFIC BACKGROUND:
--    Earth's rotation rate varies slightly due to:
--    - Tidal torques from Moon and Sun (18.6yr, 8.85yr periods)
--    - Ocean and atmospheric angular momentum exchange
--    - Core-mantle coupling
--    This function focuses on the tidal component, which is predictable
--    from celestial mechanics.
--
--  ALGORITHM:
--    1. Load dLOD time series data from file
--    2. Create forcing function from time values (linear ramp)
--    3. Convert tidal periods to frequencies (Year_Length / Period)
--    4. Perform multivariate regression to fit tidal constituents
--    5. Extract amplitude and phase for each constituent
--    6. Generate model and compute correlation coefficient
--
--  VALIDATION USE:
--    dLOD is well-understood and measurable, making it ideal for:
--    - Verifying LTE solver correctness
--    - Calibrating tidal parameters
--    - Testing regression algorithm
--    If GEM-LTE can reproduce known dLOD signals, confidence increases
--    for applying it to less-understood climate phenomena (ENSO, etc.)
--
--  RETURNS:
--    Long_Periods_Amp_Phase array containing fitted amplitude and phase
--    for each tidal constituent (18.6yr, 8.85yr, etc.)
--
--  ============================================================================

with Text_IO;
with GEM.LTE.Primitives;

function GEM.dLOD (File_Name : in String) return GEM.LTE.Long_Periods_Amp_Phase
is
   use GEM.LTE, GEM.LTE.Primitives;
   D : Data_Pairs := Make_Data (File_Name);  -- Load dLOD observations
   First, Last : Integer;
   Singular : Boolean;  -- Flags if regression matrix is singular
   Forcing : Data_Pairs := D;  -- Will hold time-based forcing function
   Model : Data_Pairs := D;    -- Will hold fitted model output
   DBLT : GEM.LTE.Long_Periods := GEM.LTE.LP; -- Tidal periods (years)
   DBLTAP : GEM.LTE.Long_Periods_Amp_Phase := GEM.LTE.LPAP;  -- Amp/Phase results
   --  TODO: Can remove - Experimental addition of annual/semi-annual harmonics
   --  was tested but standard long-period constituents proved sufficient.
   -- & Gem.LTE.Long_Periods_Amp_Phase'((0.0,0.0), (0.0,0.0), (0.0,0.0));
   Level, K0, Trend, Accel : Long_Float := 0.0;  -- Regression constants
   Last_Time : Long_Float;
begin
   First := D'First;
   Last := D'Last;
   Text_IO.Put_Line ("records from " & First'Img & " ... " & Last'Img);
   Text_IO.Put_Line
     ("factors from " & DBLT'First'Img & " ... " & DBLT'Last'Img);
   
   --  Create forcing function as linear time ramp
   --  This represents the secular trend component
   for I in Forcing'Range loop
      Forcing (I).Value := Forcing (I).Date;
      Last_Time := Forcing (I).Value;
   end loop;
   Text_IO.Put ("updated forcing  ");
   Put (Last_Time);
   Text_IO.New_Line;
   
   --  Convert periods (years) to frequencies (cycles per year)
   for I in DBLT'Range loop
      DBLT (I) := GEM.LTE.Year_Length / DBLT (I);
   end loop;

   --  Perform multivariate regression to fit tidal constituents
   Trend := 0.0;
   Regression_Factors
     (Data_Records => D, -- Time series
   --  TODO: Can remove - First/Last parameters made redundant by passing
   --  full data array. Interval selection now handled by caller if needed.
   --First => First,
   --Last => Last,  -- Training Interval

      Forcing => Forcing,  -- Value @ Time
      NM => DBLT'Last, -- # modulations
      DBLT => DBLT, --D.B.LT,
      DALTAP => DBLTAP, --D.A.LTAP,
      DALEVEL => Level, DAK0 => K0, Secular_Trend => Trend,
      Accel => Accel, Singular => Singular);

   Text_IO.Put_Line ("Singular? " & Singular'Img);
   
   --  Print fitted parameters for each constituent
   for I in DBLT'Range loop
      Put (DBLT (I));
      Text_IO.Put ("   ");
      --  NOTE: Integration step commented out - dLOD is already a rate (derivative).
      --  If analyzing absolute LOD (length of day) instead, would need:
      --    DBLTAP(I).Amplitude := DBLTAP(I).Amplitude * 26.736 / DBLT(I);
      --    DBLTAP(I).Phase := DBLTAP(I).Phase + Pi/2.0;
      Put (DBLTAP (I).Amplitude);
      Text_IO.Put ("   ");
      Put (DBLTAP (I).Phase);
      Text_IO.Put ("   ");
      Text_IO.New_Line;
   end loop;

   --  Generate model from fitted parameters and compute correlation
   Model :=
     GEM.LTE.Primitives.LTE
       (Forcing => Forcing, Wave_Numbers => DBLT, Amp_Phase => DBLTAP,
        Offset => Level, K0 => K0, Trend => 0.0, NonLin => 1.0);
   Put (CC (D, Model), "=CC ");
   Put (Year_Length, "=Yr", True);

   return DBLTAP;

end GEM.dLOD;
