--  ENSO Optimization Driver
--
--  Main procedure for optimizing Laplace's Tidal Equation (LTE) model parameters
--  to fit climate indices (e.g., ENSO/Niño data) using parallel search tasks.
--
--  The program:
--  1. Loads initial parameters and dLOD (day length) reference data
--  2. Spawns parallel worker tasks to search parameter space
--  3. Monitors progress and handles user interrupts (q/s/x keys, timeout)
--  4. Saves best-fit parameters when optimization completes

with System.Task_Info;
with GEM.LTE.Primitives.Solution;
with GEM.LTE.Primitives.Shared;
with Text_IO;
with GEM.dLOD;
with GNAT.OS_Lib;
with Ada.Calendar;
with Interfaces.C;

procedure ENSO_Opt
is  -- gprbuild lte.gpr enso_opt -largs -Wl, --stack=40000000

   --  Configuration from environment variables
   N : Positive :=  --  Number of parallel worker threads
     GEM.Getenv
       ("NUMBER_OF_PROCESSORS", System.Task_Info.Number_Of_Processors);
   dLOD_dat : String := GEM.Getenv ("DLOD_DAT", "../dlod3.dat");
   Ch : Character;  --  User keyboard input
   Avail : Boolean;  --  Keyboard input available flag
   T : Ada.Calendar.Time;  --  Timeout tracking
   Cycle : Duration :=  --  Maximum optimization duration
     Duration (GEM.Getenv ("TIMEOUT", Long_Float (Duration'Last) / 2.0));
   dLOD_Scale : constant Long_Float := GEM.Getenv ("DLOD_SCALE", 0.0);
   Expect : constant Boolean := GEM.Getenv ("EXPECT", False);
   use type Ada.Calendar.Time;

   --  D = Complete model parameter set (tidal periods, amplitudes, phases,
   --  LTE modulation parameters, impulse response coefficients)
   D : GEM.LTE.Primitives.Shared.Param_S :=
     (NLP => GEM.LTE.LP'Length, NLT => GEM.LTE.LTM'Length,
      A =>
        (k0 => 0.0, level => 0.0, NLP => GEM.LTE.LP'Length,
         NLT => GEM.LTE.LTM'Length, LP => GEM.LTE.LP, LTAP => GEM.LTE.LTAP),
      B =>
        (NLP => GEM.LTE.LP'Length, NLT => GEM.LTE.LTM'Length,
         LPAP => GEM.LTE.LPAP, LT => GEM.LTE.LTM, Offset => 0.0, bg => 0.0,
         ImpA => 1.0, ImpB => 0.0, ImpC => 0.0, DelA => 0.0, DelB => 0.0,
         Asym => 0.0, Ann1 => 0.0, Ann2 => 0.0, Sem1 => 0.0, Sem2 => 0.0,
         Year => 0.0, IR => 0.0, mA => 0.0, mP => 0.0, shiftT => 0.000_00,
         init => 0.006_3),
      C => (others => 0));

begin
   --  Stack size limit is now raised by the launcher (ulimit -s unlimited)
   --  before this executable is invoked, so the setrlimit block below is
   --  no longer needed here.
   --
   --  --  Remove the main thread stack size limit (equivalent to ulimit -s unlimited).
   --  --  The default ~8 MB limit causes STORAGE_ERROR with large data sets.
   --  declare
   --     RLIMIT_STACK : constant Interfaces.C.int := 3;
   --     type rlim_t is mod 2 ** 64;
   --     type Rlimit is record
   --        rlim_cur, rlim_max : rlim_t;
   --     end record;
   --     pragma Convention (C, Rlimit);
   --     function setrlimit
   --       (resource : Interfaces.C.int; rlp : access Rlimit)
   --        return Interfaces.C.int;
   --     pragma Import (C, setrlimit, "setrlimit");
   --     R : aliased Rlimit := (rlim_t'Last, rlim_t'Last);
   --     use type Interfaces.C.int;
   --  begin
   --     if setrlimit (RLIMIT_STACK, R'Access) /= 0 then
   --        Text_IO.Put_Line ("Warning: could not remove stack size limit");
   --     end if;
   --  end;

   declare
      --  AP = Amplitude/Phase data from day-length-of-day (dLOD) measurements
      --  Used as reference forcing for the tidal model
      AP : GEM.LTE.Long_Periods_Amp_Phase := GEM.dLOD (dLOD_dat);
   begin

      Text_IO.Put_Line (N'Img & " processors available, timeout=" & Cycle'Img);

      --  Load previously saved parameters if available (warm start)
      GEM.LTE.Primitives.Shared.Load (D);

      --  COMMENTED CODE: Year adjustment feature - disabled
      --  Would adjust the annual period parameter based on LP values
      --  TODO: Remove after confirming not needed - was experimental feature
      --  for testing sensitivity to yearly modulation
      --if GEM.Getenv("YTRIM", FALSE) then
--   GEM.LTE.Year_Adjustment(D.B.Year, D.A.LP); -- should be a protected call?
      --end if;

      Text_IO.Put_Line ("YA=" & D.B.Year'Img);

      --  Initialize tidal forcing amplitudes and phases from dLOD data
      if GEM.Command_Line_Option_Exists ("r") or
        not GEM.Getenv ("DLOD_REF", True)
      then
         --  Fresh start: use dLOD data directly with optional scaling
         Text_IO.Put_Line ("Loading dLOD");
         for I in D.B.LPAP'Range loop
            GEM.LTE.LPRef (I).Amplitude := AP (I).Amplitude;
            GEM.LTE.LPRef (I).Phase := AP (I).Phase;
            D.B.LPAP (I).Amplitude :=
              AP (I).Amplitude * (1.0 + dLOD_Scale * D.A.LP (I));
            D.B.LPAP (I).Phase := AP (I).Phase;
         end loop;
      else
         --  Reference mode: store dLOD but don't overwrite current parameters
         Text_IO.Put_Line ("Referencing dLOD");
         for I in D.B.LPAP'Range loop
            GEM.LTE.LPRef (I).Amplitude := AP (I).Amplitude;
            GEM.LTE.LPRef (I).Phase := AP (I).Phase;
         end loop;
      end if;

      --  Store initialized parameters for worker tasks to access
      GEM.LTE.Primitives.Shared.Put (D);

   end;

   T := Ada.Calendar.Clock;

   --  Launch N parallel worker tasks to search parameter space
   --  Each task performs random descent optimization independently
   GEM.LTE.Primitives.Solution.Start (D.NLP, D.NLT, N);

   --  Main monitoring loop: check status, handle user input, enforce timeout
   for I in 1 .. Integer'Last loop

      --  Get status from worker tasks (blocking call)
      declare
         S : String := GEM.LTE.Primitives.Solution.Status;
      begin
         --  Periodic progress reporting
         if I mod GEM.LTE.Primitives.Solution.Check_Every_N_Loops = 0 then
            Text_IO.Put_Line (S & " #" & I'Img);
         end if;
         exit when GEM.LTE.Primitives.Halted;
      end;

      --  Check for user keyboard commands
      if Expect then
         Text_IO.Get (Ch);  --  Blocking for scripted use
         Avail := True;
      else
         Text_IO.Get_Immediate (Ch, Avail);  --  Non-blocking for interactive
      end if;

      if Avail then
         if Ch = 'q' or Ch = 's' then
            --  'q' or 's': Stop optimization and save results
            GEM.LTE.Primitives.Stop;
         elsif Ch = 'x' then
            --  'x': Exit immediately without saving
            Text_IO.Put_Line ("Exiting, no save");
            GNAT.OS_Lib.OS_Exit (0);
         elsif Ch in '1' .. '9' then
            --  '1'-'9': Adjust correlation threshold trigger dynamically
            GEM.LTE.Primitives.Solution.Set_Trigger
              (Character'Pos (Ch) - Character'Pos ('1') + 1);
         end if;
      end if;

      if Ch = 'q' then
         exit when GEM.LTE.Primitives.Halted;
      end if;

      --  Enforce timeout limit
      if Ada.Calendar.Clock > T + Cycle then
         GEM.LTE.Primitives.Stop;
         T := Ada.Calendar.Clock;
         exit when GEM.LTE.Primitives.Halted;
      end if;

   end loop;

   --  COMMENTED CODE: Debug message for main thread exit
   --  TODO: Can remove - leftover from debugging parallel task shutdown
   --Text_IO.Put_Line("Main exiting, flushing other tasks");

   --  Wait for worker tasks to complete their final output writes
   delay 5.0;

end ENSO_Opt;
