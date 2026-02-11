with Text_IO;
with GEM.LTE.Primitives;

function GEM.Mix_Regression
  (Target : in LTE.Primitives.Data_Pairs; Start, Finish : in Long_Float;
   First_Index, Last_Index : in Integer := 0; Exclude : in Boolean := False;
   Display_Results : in Boolean := False) return LTE.Primitives.Data_Pairs
is
   use GEM.LTE, GEM.LTE.Primitives;
   FT : Text_IO.File_Type;
   D : Data_Pairs := Target; -- := Make_Data(File_Name);
   First, Last : Integer;
   Singular : Boolean;
   Forcing : Data_Pairs := D;
   Model : Data_Pairs := D;
   P_A1, P_A2, P_A3, P_A4, P_A5, P_A6, P_A7, P_A8, P_A9, P_A10, P_A11,
   P_A12 : GEM.LTE.Long_Periods := GEM.LTE.LP_Set;
   P_B1 : GEM.LTE.Long_Periods := GEM.LTE.LP_RSet;
   P_Q1 : GEM.LTE.Long_Periods := GEM.LTE.LP_QSet;
   P_Annual : GEM.LTE.Long_Periods := GEM.LTE.LP_Annual;
   A_A1, A_A2, A_A3, A_A4, A_A5, A_A6, A_A7, A_A8, A_A9, A_A10, A_A11,
   A_A12 : GEM.LTE.Long_Periods_Amp_Phase := GEM.LTE.LPAP_Set;
   A_B1 : GEM.LTE.Long_Periods_Amp_Phase := GEM.LTE.LPAP_RSet;
   A_Q1 : GEM.LTE.Long_Periods_Amp_Phase := GEM.LTE.LPAP_QSet;
   A_Annual : GEM.LTE.Long_Periods_Amp_Phase := GEM.LTE.LPAP_Annual;
   Level, K0 : Long_Float;
   Trend : Long_Float := 1.0;
   Accel : Long_Float;
   Last_Time : Long_Float;
   Freq, Freq2 : Long_Float;
   One : constant Long_Float := 1.0;
   function Match (Freq : in Long_Float) return Boolean is
      M : constant Long_Float := 6.0;
   begin
      for I in P_A1'Range loop -- to aliased frequency
         if Integer (M * P_A1 (I)) = Integer (M * Freq) then
            if Freq < 0.07 then
               return True;
            else
               return False;
            end if;
         end if;
      end loop;
      for I in P_B1'Range loop -- to aliased frequency
         if Integer (M * P_B1 (I)) = Integer (M * Freq) then
            if Freq < 0.07 then
               return True;
            else
               return False;
            end if;
         end if;
      end loop;
      return False;
   end Match;

   function Index (Fraction : in Long_Float) return Positive is
   begin
      return First + abs Integer (Fraction * Long_Float (Last - First));
   end Index;

   function Interval
     (TS : in Data_Pairs; Complement : in Boolean := False) return Data_Pairs
   is
      S, F : Integer;
   begin
      if First_Index /= 0 then
         S := abs First_Index;
      else
         S := First + abs Integer (Start * Long_Float (Last - First));
      end if;
      if Last_Index /= 0 then
         F := abs Last_Index;
      else
         F := First + abs Integer (Finish * Long_Float (Last - First));
      end if;
      if not Complement then
         if First_Index > 0 then
            return TS (S .. F);
         elsif First_Index < 0 then
            return TS (TS'First .. S) & TS (F .. TS'Last);
         elsif Start > 0.0 then
            --Text_IO.Put_Line("Fit from " & S'Img & " .. " & F'Img);
            return TS (S .. F);
         else
            --Text_IO.Put_Line("Fit from 1 .. " & S'Img & " and " & F'Img & " .. " & Last'Img);
            return TS (TS'First .. S) & TS (F .. TS'Last);
         end if;
      else
         if First_Index > 0 then
            return TS (TS'First .. S) & TS (F .. TS'Last);
         elsif First_Index > 0 then
            return TS (S .. F);
         elsif Start > 0.0 then
            --Text_IO.Put_Line("Exclude fit from 1 .. " & S'Img & " and " & F'Img & " .. " & Last'Img);
            return TS (TS'First .. S) & TS (F .. TS'Last);
         else
            --Text_IO.Put_Line("Exclude Fit from " & S'Img & " .. " & F'Img);
            return TS (S .. F);
         end if;
      end if;
   end Interval;

   procedure Fit
     (P : in GEM.LTE.Long_Periods; A : in out GEM.LTE.Long_Periods_Amp_Phase)
   is
   --             Start, Finish : in Long_Float;
   --             Complement : in Boolean := False) is

   begin
      Regression_Factors
        (Data_Records => Interval (D), -- Time series
         Forcing => Interval (Forcing),  -- Value @ Time
         NM => P'Last, -- # modulations
         DBLT => P, --D.B.LT,
         DALTAP => A, --D.A.LTAP,
         DALEVEL => Level, DAK0 => K0,
         Secular_Trend => Trend, Accel => Accel, Singular => Singular);
   end Fit;

   function Ex (A : in GEM.LTE.Long_Periods) return GEM.LTE.Long_Periods is
   begin
      if Exclude then
         return A (0 .. -1);
      else
         return A;
      end if;
   end Ex;

   function Ex
     (A : in GEM.LTE.Long_Periods_Amp_Phase)
      return GEM.LTE.Long_Periods_Amp_Phase
   is
   begin
      if Exclude then
         return A (0 .. -1);
      else
         return A;
      end if;
   end Ex;

begin
   First := D'First;
   Last := D'Last;
--   Text_IO.Put_Line("records from " & First'Img & " ... " & Last'Img);
--   Text_IO.Put_Line("factors from " & P_A1'First'Img & " ... " & P_A1'Last'Img);
   for I in Forcing'Range loop
      Forcing (I).Value := Forcing (I).Date;
      Last_Time := Forcing (I).Value;
   end loop;
--   Text_IO.Put("updated forcing  ");
--   Put(Last_Time);
--   Text_IO.New_Line;
   for I in P_A1'Range loop -- to aliased frequency
      Freq := GEM.LTE.Year_Length / P_A1 (I);
      Freq := abs Long_Float'Remainder (Freq, One);
      P_A1 (I) := Freq;
      P_A2 (I) := One - Freq;
      P_A3 (I) := One + Freq;
      P_A4 (I) := 2.0 * One - Freq;
      P_A5 (I) := 2.0 * One + Freq;
      P_A6 (I) := 3.0 * One - Freq;
      P_A7 (I) := 4.0 * One - Freq;
      P_A8 (I) := 5.0 * One - Freq;
      P_A9 (I) := 3.0 * One + Freq;
      P_A10 (I) := 4.0 * One + Freq;
      P_A11 (I) := 5.0 * One + Freq;
      P_A12 (I) := 6.0 * One - Freq;
   end loop;

   for I in P_B1'Range loop -- to aliased frequency
      Freq := 2.0 * GEM.LTE.Year_Length / P_B1 (I);

      if Integer (Freq) mod 2 = 1 then
         Freq := abs Long_Float'Remainder (Freq, One);
      else
         Freq := One - abs Long_Float'Remainder (Freq, One);
      end if;
      P_B1 (I) := Freq * 0.5;

--      if Match(Freq) then
      -- P_B1(I) := 100.0*(Long_Float(200*I) + 313.131313/Long_Float(I));
--         P_B1(I) := 3.0 + Freq;
--      else
--         P_B1(I) := Freq;
--      end if;
--      if Match(Freq2) then
--         -- P_B2(I) := 1000.0*(Long_Float(100*I)  + 313.131313/Long_Float(I));
--         P_B2(I) := 3.0 + Freq2;
--      else
--         P_B2(I) := Freq2;
--      end if;

      --P_B1(I) := Freq*0.5;
      -- P_B2(I) := (One-Freq)*0.5;

--      if Integer(Freq2) mod 2 = 1 then
--         Freq2 := abs Long_Float'Remainder(Freq2,One);
--      else
--         Freq2 := One - abs Long_Float'Remainder(Freq2,One);
--      end if;
--      P_Q1(I) := Freq2*0.25;
   end loop;

   for I in P_Q1'Range loop -- to aliased frequency
      Freq2 := 4.0 * GEM.LTE.Year_Length / P_Q1 (I);
      if I = P_Q1'First then
         Freq2 := 0.25 * (One - abs Long_Float'Remainder (Freq2, One));
      else
         Freq2 := 0.25 * (abs Long_Float'Remainder (Freq2, One));
      end if;

--      if Match(Freq2) then
--         P_Q1(I) := 4.0 - Freq2;
--      else
      P_Q1 (I) := Freq2;
--      end if;
   end loop;

   declare
      P : GEM.LTE.Long_Periods :=
        P_A1 &
        P_A2 &
        P_A3 &
        P_A4 &
        P_A5 &
        P_A6 &
        P_A7 &
        P_A8   -- & P_A9 & P_A10 & P_A11 & P_A12
      &
        Ex (P_B1 & P_Q1) & P_Annual; -- & P_Q1;
      A : GEM.LTE.Long_Periods_Amp_Phase :=
        A_A1 &
        A_A2 &
        A_A3 &
        A_A4 &
        A_A5 &
        A_A6 &
        A_A7 &
        A_A8  --   & A_A9 & A_A10 & A_A11 & A_A12
      &
        Ex (A_B1 & A_Q1) & A_Annual; -- & A_Q1;

      procedure Print is
      begin
         if Display_Results then
            Text_IO.Put_Line ("Singular? " & Singular'Img);
            for I in P'Range loop
               Put (P (I));
               Text_IO.Put ("   ");
               -- INTEGRATE dLOD -- DBLTAP(I).Amplitude := DBLTAP(I).Amplitude * 26.736 / DBLT(I);
               Put (A (I).Amplitude);
               Text_IO.Put ("   ");
-- INTEGRATE dLOD -- DBLTAP(I).Phase := DBLTAP(I).Phase + Ada.Numerics.Pi/2.0;
               Put (A (I).Phase);
               Text_IO.Put ("   ");
               Text_IO.New_Line;
            end loop;
         end if;
      end Print;
      procedure Debug is
      begin
         if Display_Results then
            Put (CC (D, Model), "=CC ");
            Put (CC (Interval (D), Interval (Model)), "=CC(training) ");
            Put
              (CC (Interval (D, True), Interval (Model, True)),
               "=CC(validate) ");
            Put
              (DTW_Distance (Interval (D, True), Interval (Model, True), 1),
               "=DTW(validate) ");

            Put (Year_Length, "=Yr ", True);
            Put (Model (Index (Start)).Date, "=start ");
            Put (Model (Index (Finish)).Date, "=finish ");

            GEM.LTE.Primitives.Save (Model, D, Forcing);

            Text_IO.Create (FT, Text_IO.Out_File, "training.txt");
            Text_IO.Put_Line (FT, Model (Index (Start)).Date'Img & ", 0.0");
            Text_IO.Put_Line (FT, Model (Index (Finish)).Date'Img & ", 0.0");
            Text_IO.Close (FT);
         end if;
      end Debug;

   begin
      Trend := 0.0;
      Fit (P, A); --, Start, Finish);

      Print;

      Model :=
        GEM.LTE.Primitives.LTE
          (Forcing => Forcing, Wave_Numbers => P, Amp_Phase => A,
           Offset => Level, K0 => K0, Trend => Trend, Accel => Accel,
           NonLin => 1.0);

      Debug;
      return Model;
      --  declare
      --     PS : constant Gem.LTE.Period_Set := (N => P'Length,
      --                              LP => P,
      --                              AP => A);
      --  begin
      --     return PS;
      --  end;
      --return A;
   end;

end GEM.Mix_Regression;
