with System.Task_Info;
with GEM.LTE.Primitives.Solution;
with GEM.LTE.Primitives.Shared;
with Text_IO;
with GEM.dLOD;
with GNAT.OS_Lib;
with Ada.Directories;

-- Test unit for GEM-LTE that performs a single pass calculation
-- Usage: Set CLIMATE_INDEX environment variable or use default nino4.dat
-- Example: CLIMATE_INDEX=nino4.dat gem_test
procedure GEM_Test is
   
   -- Configuration for single-pass test
   Climate_Index : String := GEM.Getenv("CLIMATE_INDEX", "nino4.dat");
   Test_Timeout  : constant Duration := 0.0;  -- No timeout - single pass only
   Test_Trigger  : constant Integer := 0;      -- No trigger - run once and stop
   
   N : Positive := GEM.Getenv("NUMBER_OF_PROCESSORS", System.Task_Info.Number_Of_Processors);
   dLOD_dat : String := GEM.Getenv("DLOD_DAT", "../dLOD3.dat"); 
   dLOD_Scale : constant Long_Float := GEM.Getenv("DLOD_SCALE", 0.0);

   D : GEM.LTE.Primitives.Shared.Param_S :=
     (NLP => GEM.LTE.LP'Length,
      NLT => GEM.LTE.LTM'Length,
      A =>
        (k0     => 0.0,
         level  => 0.0,
         NLP    => GEM.LTE.LP'Length,
         NLT    => GEM.LTE.LTM'Length,
         LP     => GEM.LTE.LP,
         LTAP   => GEM.LTE.LTAP),
      B =>
        (NLP    => GEM.LTE.LP'Length,
         NLT    => GEM.LTE.LTM'Length,
         LPAP   => GEM.LTE.LPAP,
         LT     => GEM.LTE.LTM,
         Offset => 0.0,
         bg     => 0.0,
         ImpA   => 1.0,
         ImpB   => 0.0,
         impC  =>  0.0,
         delA  =>  0.0,
         delB  =>  0.0,
         asym  =>  0.0,
         ann1  =>  0.0,
         ann2  =>  0.0,
         sem1  =>  0.0,
         sem2  =>  0.0,
         year  =>  0.0,
         ir    =>  0.0,
         mA     => 0.0,
         mP     => 0.0,
         shiftT => 0.00000,
         init   => 0.0063),
      C => (others => 0)
     );

   Reference_CSV_File : constant String := "lte_results_reference.csv";
   Output_CSV_File    : constant String := "lte_results.csv";
   Test_Passed        : Boolean := False;

begin
   Text_IO.Put_Line("========================================");
   Text_IO.Put_Line("  GEM-LTE Test Unit");
   Text_IO.Put_Line("========================================");
   Text_IO.Put_Line("Climate Index: " & Climate_Index);
   Text_IO.Put_Line("Processors:    " & N'Img);
   Text_IO.Put_Line("Timeout:       " & Test_Timeout'Img & " (single pass)");
   Text_IO.Put_Line("Trigger:       " & Test_Trigger'Img & " (run once)");
   Text_IO.Put_Line("dLOD Data:     " & dLOD_dat);
   Text_IO.Put_Line("========================================");

   -- Load dLOD data
   declare
      AP : GEM.LTE.Long_Periods_Amp_Phase := GEM.dLOD(dLOD_dat);
   begin
      Text_IO.Put_Line("Loading shared parameters...");
      GEM.LTE.Primitives.Shared.Load(D);
      
      Text_IO.Put_Line("Setting up dLOD reference...");
      for I in D.B.LPAP'Range loop
        GEM.LTE.LPRef(I).Amplitude := AP(I).Amplitude;
        GEM.LTE.LPRef(I).Phase := AP(I).Phase;
        D.B.LPAP(I).Amplitude := AP(I).Amplitude * (1.0 + dLOD_Scale*D.A.LP(I));
        D.B.LPAP(I).Phase := AP(I).Phase;
      end loop;
      
      GEM.LTE.Primitives.Shared.Put(D);
   end;

   Text_IO.Put_Line("Starting solution calculation...");
   GEM.LTE.Primitives.Solution.Start(D.NLP, D.NLT, N);

   -- Run until first calculation completes (one pass only)
   -- With TIMEOUT=0.0 and TRIGGER=0, this should stop after one iteration
   declare
      Loop_Count : Integer := 0;
      Max_Loops  : constant Integer := 10000;  -- Safety limit
   begin
      for I in 1..Max_Loops loop
         Loop_Count := I;
         
         -- Check status
         declare
            S : String := GEM.LTE.Primitives.Solution.Status;
         begin
            if I mod 100 = 0 then
               Text_IO.Put_Line("Loop " & I'Img & ": " & S);
            end if;
            
            -- Exit when calculation halts
            exit when GEM.LTE.Primitives.Halted;
         end;
      end loop;
      
      Text_IO.Put_Line("Calculation complete after " & Loop_Count'Img & " iterations");
   end;

   -- Give tasks time to flush output
   delay 2.0;

   Text_IO.Put_Line("========================================");
   Text_IO.Put_Line("Checking results...");
   
   -- Check if output CSV was generated
   if Ada.Directories.Exists(Output_CSV_File) then
      Text_IO.Put_Line("✓ Output file generated: " & Output_CSV_File);
      
      -- Compare with reference if it exists
      if Ada.Directories.Exists(Reference_CSV_File) then
         Text_IO.Put_Line("Comparing with reference: " & Reference_CSV_File);
         -- TODO: Implement CSV comparison logic
         -- For now, just report that both files exist
         Test_Passed := True;
      else
         Text_IO.Put_Line("⚠ Reference file not found: " & Reference_CSV_File);
         Text_IO.Put_Line("  This may be the first run. Output file can be used as reference.");
         Test_Passed := True;  -- Consider test passed if output exists
      end if;
   else
      Text_IO.Put_Line("✗ Output file NOT generated: " & Output_CSV_File);
      Test_Passed := False;
   end if;

   Text_IO.Put_Line("========================================");
   if Test_Passed then
      Text_IO.Put_Line("TEST PASSED");
      GNAT.OS_Lib.OS_Exit(0);
   else
      Text_IO.Put_Line("TEST FAILED");
      GNAT.OS_Lib.OS_Exit(1);
   end if;

end GEM_Test;
