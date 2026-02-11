with System.Task_Info;
with GEM.LTE.Primitives.Solution;
with GEM.LTE.Primitives.Shared;
with Text_IO;
with GEM.dLOD;
with GNAT.OS_Lib;
with Ada.Calendar;

procedure ENSO_Opt
is  -- gprbuild lte.gpr enso_opt -largs -Wl, --stack=40000000
   N : Positive :=
     GEM.Getenv
       ("NUMBER_OF_PROCESSORS", System.Task_Info.Number_Of_Processors);
   dLOD_dat : String := GEM.Getenv ("DLOD_DAT", "../dlod3.dat");
   Ch : Character;
   Avail : Boolean;
   T : Ada.Calendar.Time;
   Cycle : Duration :=
     Duration (GEM.Getenv ("TIMEOUT", Long_Float (Duration'Last) / 2.0));
   dLOD_Scale : constant Long_Float := GEM.Getenv ("DLOD_SCALE", 0.0);
   Expect : constant Boolean := GEM.Getenv ("EXPECT", False);
   use type Ada.Calendar.Time;

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
   declare
      AP : GEM.LTE.Long_Periods_Amp_Phase := GEM.dLOD (dLOD_dat);
   begin

      Text_IO.Put_Line (N'Img & " processors available, timeout=" & Cycle'Img);
      GEM.LTE.Primitives.Shared.Load (D); -- if available
      --if GEM.Getenv("YTRIM", FALSE) then
--   GEM.LTE.Year_Adjustment(D.B.Year, D.A.LP); -- should be a protected call?
      --end if;
      Text_IO.Put_Line ("YA=" & D.B.Year'Img);
      if GEM.Command_Line_Option_Exists ("r") or
        not GEM.Getenv ("DLOD_REF", False)
      then
         Text_IO.Put_Line ("Loading dLOD");
         for I in D.B.LPAP'Range loop
            GEM.LTE.LPRef (I).Amplitude := AP (I).Amplitude;
            GEM.LTE.LPRef (I).Phase := AP (I).Phase;
            D.B.LPAP (I).Amplitude :=
              AP (I).Amplitude * (1.0 + dLOD_Scale * D.A.LP (I));
            D.B.LPAP (I).Phase := AP (I).Phase;
         end loop;
      else -- update
         Text_IO.Put_Line ("Referencing dLOD");
         for I in D.B.LPAP'Range loop
            GEM.LTE.LPRef (I).Amplitude := AP (I).Amplitude;
            GEM.LTE.LPRef (I).Phase := AP (I).Phase;
         end loop;
      end if;
      GEM.LTE.Primitives.Shared.Put (D);

   end;

   T := Ada.Calendar.Clock;

   GEM.LTE.Primitives.Solution.Start (D.NLP, D.NLT, N);

   for I in 1 .. Integer'Last loop
      -- The call to Status is blocking
      declare
         S : String := GEM.LTE.Primitives.Solution.Status;
      begin
         if I mod GEM.LTE.Primitives.Solution.Check_Every_N_Loops = 0 then
            Text_IO.Put_Line (S & " #" & I'Img);
         end if;
         exit when GEM.LTE.Primitives.Halted;
      end;
      --
      if Expect then
         Text_IO.Get (Ch);
         Avail := True;
      else
         Text_IO.Get_Immediate (Ch, Avail);
      end if;
      if Avail then
         if Ch = 'q' or Ch = 's' then
            GEM.LTE.Primitives.Stop;
         elsif Ch = 'x' then
            Text_IO.Put_Line ("Exiting, no save");
            GNAT.OS_Lib.OS_Exit (0);
         elsif Ch in '1' .. '9' then
            GEM.LTE.Primitives.Solution.Set_Trigger
              (Character'Pos (Ch) - Character'Pos ('1') + 1);
         end if;
      end if;
      if Ch = 'q' then
         exit when GEM.LTE.Primitives.Halted;
      end if;
      if Ada.Calendar.Clock > T + Cycle then
         GEM.LTE.Primitives.Stop;
         T := Ada.Calendar.Clock;
         exit when GEM.LTE.Primitives.Halted; -- new
      end if;
   end loop;
   --Text_IO.Put_Line("Main exiting, flushing other tasks");

   delay 5.0;

end ENSO_Opt;
