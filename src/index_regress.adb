with Text_IO;
with GEM.LTE.Primitives;
with GEM.Mix_Regression;
with Ada.Command_Line;

procedure index_regress is
   FT : Text_IO.File_Type;
begin

   declare
      Name : constant String := Ada.Command_Line.Argument (1);
      Flts : constant GEM.Fs := GEM.S_to_LF (Ada.Command_Line.Argument (2));
      -- "dlod3.dat"
      Target : GEM.LTE.Primitives.Data_Pairs :=
        GEM.LTE.Primitives.Make_Data (Name);
      Result : GEM.LTE.Primitives.Data_Pairs :=
        GEM.Mix_Regression
          (Target, Flts (1), Flts (2), 0, 0, GEM.Getenv ("QUADEXCLUDE", False),
           True);

--      Freqs : Gem.Lte.Periods := Ps.LP;
--      TS : Gem.LTE.Primitives.Data_Pairs(1..1000) := (others => (0.0, 0.0));
--      Res : Gem.LTE.Primitives.Data_Pairs := TS;
--      T : Long_FLoat := 1900.0;
   begin
      --  for I in TS'Range loop
      --     TS(I).Date := T;
      --     T := T + 0.001;
      --  end loop;
      --  for I in Freqs'Range loop
      --     Freqs(I) := 6.28/Freqs(I);
      --  end loop;
      --  Res := Gem.LTE.primitives.Tide_Sum(TS, PS.AP, Freqs);
      --  Text_IO.Create(FT, Text_IO.Out_File,  "fine_results.csv");
      --  for I in Res'Range loop
      --     Text_IO.Put_Line(FT, Res(I).Date'Img & " " & Res(I).Value'Img);
      --  end loop;
      --  Text_IO.Close(FT);

      Text_IO.Create (FT, Text_IO.Out_File, "lte_label.txt");
      Text_IO.Put_Line (FT, Name);
      if Flts (1) < 0.0 then
         Text_IO.Put_Line (FT, "Fraction Exclude");
      else
         Text_IO.Put_Line (FT, "Fraction Include");
      end if;
      Text_IO.Put_Line (FT, Ada.Command_Line.Argument (2));
      Text_IO.Close (FT);

      -- return Result;

   end;

end index_regress;

-- 0.99064132219=CC  365.24123840000=Yr
-- 0.99064647703=CC  365.24223840000=Yr
-- 0.99064997996=CC  365.24323840000=Yr
-- 0.99065181974=CC  365.24423840000=Yr
-- 0.99065198518=CC  365.24523840000=Yr  ***
-- 0.99065046507=CC  365.24623840000=Yr
-- 0.99064724829=CC  365.24723840000=Yr
