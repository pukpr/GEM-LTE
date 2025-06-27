with Text_IO;
with GEM.LTE.Primitives;

function GEM.Mix_Regression (File_Name : in String) return Gem.LTE.Long_Periods_Amp_Phase is
   use GEM.LTE, GEM.LTE.Primitives;
   D : Data_Pairs := Make_Data(File_Name);
   First, Last : Integer;
   Singular : Boolean;
   Forcing : Data_Pairs := D;
   Model : Data_Pairs := D;
   P_A1, P_A2, P_A3, P_A4, P_A5 : Gem.LTE.Long_Periods := Gem.LTE.LP_Set;
   P_B1, P_Q1                   : Gem.LTE.Long_Periods := Gem.LTE.LP_RSet;
   A_A1, A_A2, A_A3, A_A4, A_A5 : Gem.LTE.Long_Periods_Amp_Phase := Gem.LTE.LPAP_Set; 
   A_B1, A_Q1                   : Gem.LTE.Long_Periods_Amp_Phase := Gem.LTE.LPAP_RSet;
   Level, K0, Trend, Accel : Long_Float := 0.0;
   Last_Time : Long_Float;
   Freq, Freq2 : Long_Float;
   One : constant Long_Float := 1.0;
begin
   First := D'First;
   Last := D'Last;
   Text_IO.Put_Line("records from " & First'Img & " ... " & Last'Img);
   Text_IO.Put_Line("factors from " & P_A1'First'Img & " ... " & P_A1'Last'Img);
   for I in Forcing'Range loop
      Forcing(I).Value := Forcing(I).Date;
      Last_Time := Forcing(I).Value;
   end loop;
   Text_IO.Put("updated forcing  ");
   Put(Last_Time);
   Text_IO.New_Line;
   for I in P_A1'Range loop -- to aliased frequency
      Freq := GEM.LTE.Year_Length/P_A1(I);
      Freq := abs Long_Float'Remainder(Freq,One);
      P_A1(I) := Freq;
      P_A2(I) := One - Freq;
      P_A3(I) := One + Freq;
      P_A4(I) := 2.0*One - Freq;
      P_A5(I) := 2.0*One + Freq;
   end loop;

   for I in P_B1'Range loop -- to aliased frequency
      Freq := 2.0*GEM.LTE.Year_Length/P_B1(I);
      Freq2 := 4.0*GEM.LTE.Year_Length/P_B1(I);
      Text_IO.Put_Line("Q=" & Freq2'Image);
      --Text_IO.Put_Line("F=" & Freq'Image);
      if Integer(Freq) mod 2 = 1 then
         Freq := abs Long_Float'Remainder(Freq,One);
      else
         Freq := One - abs Long_Float'Remainder(Freq,One);
      end if;
      P_B1(I) := Freq*0.5;

      if Integer(Freq2) mod 2 = 1 then
         Freq2 := abs Long_Float'Remainder(Freq2,One);
      else
         Freq2 := One - abs Long_Float'Remainder(Freq2,One);
      end if;
      P_Q1(I) := Freq2*0.25;
   end loop;

   declare
      P : Gem.LTE.Long_Periods := P_A1 & P_A2 & P_A3 & P_A4 & P_A5 & P_B1; -- & P_Q1;
      A : Gem.LTE.Long_Periods_Amp_Phase := A_A1 & A_A2 & A_A3 & A_A4 & A_A5 & A_B1; -- & A_Q1;
   begin
   Trend := 0.0;
   Regression_Factors (Data_Records => D(1..800), -- Time series
                       Forcing => Forcing(1..800),  -- Value @ Time
                       NM => P'Last, -- # modulations
                       DBLT => P, --D.B.LT,
                       DALTAP => A, --D.A.LTAP,
                       DALEVEL => LEVEL,
                       DAK0 => K0,
                       Secular_Trend => Trend,
                       Accel => Accel,
                       Singular => Singular
                       );

   Text_IO.Put_Line("Singular? " & Singular'Img);
   for I in P'Range loop
      Put(P(I)); Text_IO.Put("   ");
      -- INTEGRATE dLOD -- DBLTAP(I).Amplitude := DBLTAP(I).Amplitude * 26.736 / DBLT(I);
      Put(A(I).Amplitude); Text_IO.Put("   ");
      -- INTEGRATE dLOD -- DBLTAP(I).Phase := DBLTAP(I).Phase + Ada.Numerics.Pi/2.0;
      Put(A(I).Phase); Text_IO.Put("   ");
      Text_IO.New_Line;
   end loop;

   Model := Gem.LTE.Primitives.LTE(Forcing => Forcing,
                Wave_Numbers => P,
                Amp_Phase => A,
                Offset => Level,
                K0 => K0,
                Trend => 0.0,
                NonLin => 1.0);
   Put(CC(D,Model), "=CC " );
      Put(Year_Length, "=Yr", True);
      
   GEM.LTE.Primitives.Save(Model, D, Forcing);

   return A;
   end;

end;
