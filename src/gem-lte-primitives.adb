with Ada.Long_Float_Text_IO;
with Text_IO;
with Ada.Numerics.Long_Elementary_Functions;
with Ada.Exceptions;
with Ada.Numerics.Generic_Real_Arrays;
with GNAT.OS_Lib;
--with GEM.Matrices;

package body GEM.LTE.Primitives is
   Aliased_Period : constant Boolean := GEM.Getenv ("ALIAS", False);
   Min_Entropy : constant Boolean := GEM.Getenv ("METRIC", "") = "ME";
   Linear_Step : constant Boolean := GEM.Getenv ("STEP", False);
   Every_N_Line : constant Integer := GEM.Getenv ("EVERY", 1);
   Magnify : constant Integer := GEM.Getenv ("MAG", 1);
   Clip : constant Long_Float := GEM.Getenv ("CLIP", 0.0);
   Start_Year : constant Long_Float := GEM.Getenv ("CC_START", 0.0);
   End_Year : constant Long_Float := GEM.Getenv ("CC_END", 999_999_999.0);
   Sinc : constant Long_Float := GEM.Getenv ("SINC", 0.0);
   Modulation : constant Boolean := GEM.Getenv ("MODULATION", True); -- FALSE
   Ealign : constant Integer := GEM.Getenv ("EALIGN", 24);
   Custom_Tide : constant Boolean := GEM.Getenv ("CUSTOM", False);

   function Is_Minimum_Entropy return Boolean is
   begin
      return Min_Entropy;
   end Is_Minimum_Entropy;

   function Reduce
     (Raw : in Data_Pairs; -- Raw starts with line 1
      Every : in Positive) return Data_Pairs
   is
      Res : Data_Pairs (1 .. Raw'Length / Every) := (others => (0.0, 0.0));
   begin
      for I in Res'Range loop
         for J in I * Every .. I * Every + Every - 1 loop
            Res (I).Value := Res (I).Value + Raw (J - Every + 1).Value;
         end loop;
         Res (I).Value := Res (I).Value / Long_Float (Every);
         Res (I).Date := Raw (I * Every - Every / 2).Date;
      end loop;
      return Res;
   end Reduce;

   function Expand
     (Raw : in Data_Pairs; -- Raw starts with line 1
      Mag : in Positive) return Data_Pairs
   is
      Res : Data_Pairs (1 .. Raw'Length * Mag) := (others => (0.0, 0.0));
      Stride : Integer := 1;
   begin
      for I in Raw'First .. Raw'Last - 1 loop
         for J in 0 .. Mag - 1 loop
            Res (Stride + J).Value :=
              Raw (I).Value +
              (Raw (I + 1).Value - Raw (I).Value) * Long_Float (J) /
                Long_Float (Mag);
            Res (Stride + J).Date :=
              Raw (I).Date +
              (Raw (I + 1).Date - Raw (I).Date) * Long_Float (J) /
                Long_Float (Mag);
         end loop;
         Stride := Stride + Mag;
      end loop;
      Res (Raw'Length * Mag).Value := Raw (Raw'Last).Value;
      Res (Raw'Length * Mag).Date := Raw (Raw'Last).Date;
      return Res;
   end Expand;

   function File_Lines (Name : String) return Integer is
      Data : Text_IO.File_Type;
      Count : Integer := 0;
   begin
      if Name = "" then
         return 0;
      end if;
      Text_IO.Open (File => Data, Mode => Text_IO.In_File, Name => Name);
      loop
         Text_IO.Skip_Line (Data);
         Count := Count + 1;
      end loop;
   exception
      when Text_IO.End_Error =>
         Text_IO.Close (File => Data);
         return Count;
      when Text_IO.Name_Error =>
         Text_IO.Put_Line (Name & " not found?");
         raise Text_IO.Name_Error;
      when E : others =>
         Text_IO.Put_Line (Ada.Exceptions.Exception_Information (E));
         return 0;
   end File_Lines;

   function Make_Data (Name : String) return Data_Pairs is
      Lines : Integer := File_Lines (Name);
      Data : Text_IO.File_Type;
      Date, Value : Long_Float;
      Arr : Data_Pairs (1 .. Lines);
      Val : Long_Float;
   begin
      if Lines = 0 then
         return Arr;
      end if;
      Text_IO.Open (File => Data, Mode => Text_IO.In_File, Name => Name);
      for I in 1 .. Lines loop
         Ada.Long_Float_Text_IO.Get (File => Data, Item => Date);
         Ada.Long_Float_Text_IO.Get (File => Data, Item => Value);
         Text_IO.Skip_Line (Data);
         -- Text_IO.Put_Line(Date'Img & " " & Value'Img);
         if Clip > 0.0 then
            if Value > Clip then
               Val := Clip;
            elsif Value < -Clip then
               Val := -Clip;
            else
               Val := Value;
            end if;
            Arr (I) := (Date, Val);
         else
            Arr (I) := (Date, Value);
         end if;
      end loop;
      Text_IO.Close (File => Data);
      if Every_N_Line = 1 then
         return Arr;
      else
         return Reduce (Arr, Every_N_Line);
      end if;
   exception
      when E : others =>
         Text_IO.Put_Line (Ada.Exceptions.Exception_Information (E));
         GNAT.OS_Lib.OS_Exit (0);
         return Arr;
   end Make_Data;

   -- a running mean filter, used to damp the inpulse response
   function IIR
     (Raw : in Data_Pairs; lagA, lagB, lagC : in Long_Float;
      iA, iB, iC : in Long_Float := 0.0;
      Start : in Long_Float := Long_Float'First; mA, mB : in Long_Float := 0.0)
      return Data_Pairs
   is
      Res : Data_Pairs := Raw;
      Start_Index : Integer;
      Ramp : Long_Float;

   begin
      for I in Raw'Range loop
         Start_Index := I;
         exit when Raw (I).Date > Start;
      end loop;

      -- The first few values are sensitive to prior information so can adjust
      -- these for better fits in the early part of the time series.
      Res (Start_Index).Value := iA;
      for I in Start_Index + 1 .. Raw'Last loop
         Ramp := Long_Float'Copy_Sign (lagC, Res (I - 1).Value);
         Res (I).Value := Raw (I).Value + Res (I - 1).Value - Ramp;
      end loop;

   -- If Start is inside the series, create "pre-history" by running backwards
      for I in reverse (Raw'First + 1) .. Start_Index loop
         Ramp := Long_Float'Copy_Sign (lagC, Res (I).Value);
         Res (I - 1).Value := -Raw (I - 1).Value + Res (I).Value + Ramp;
      end loop;

      return Res;
   end IIR;

   -- a running mean filter, used to damp the inpulse response
   function IIR_1
     (Raw : in Data_Pairs; lagA, lagB, lagC : in Long_Float;
      iA, iB, iC : in Long_Float := 0.0;
      Start : in Long_Float := Long_Float'First; mA, mB : in Long_Float := 0.0)
      return Data_Pairs
   is
      Res : Data_Pairs := Raw;
      Start_Index : Integer;
      Ramp : Long_Float;

   begin
      for I in Raw'Range loop
         Start_Index := I;
         exit when Raw (I).Date > Start;
      end loop;

      -- The first few values are sensitive to prior information so can adjust
      -- these for better fits in the early part of the time series.
      Res (Start_Index).Value := 0.0;
      Res (Start_Index + 1).Value := 0.0;
      Res (Start_Index + 2).Value := iA;
      Res (Start_Index + 3).Value := iA;
      for I in Start_Index + 4 .. Raw'Last loop
         Ramp := Long_Float'Copy_Sign (lagC, Res (I - 1).Value);
         Res (I).Value := Raw (I).Value + Res (I - 1).Value - Ramp;
      end loop;

   -- If Start is inside the series, create "pre-history" by running backwards
      for I in reverse (Raw'First + 2) .. (Start_Index + 3) loop
         Ramp := Long_Float'Copy_Sign (lagC, Res (I).Value);
         Res (I - 1).Value := -Raw (I - 1).Value + Res (I).Value + Ramp;
      end loop;

      return Res;
   end IIR_1;

   function IIR_0
     (Raw : in Data_Pairs; lagA, lagB, lagC : in Long_Float;
      iA, iB, iC : in Long_Float := 0.0;
      Start : in Long_Float := Long_Float'First; mA, mB : in Long_Float := 0.0)
      return Data_Pairs
   is
      Res : Data_Pairs := Raw;
      Start_Index : Integer;
      Ramp : Long_Float;
      FB2 : Long_Float;
      Pi : Long_Float := Ada.Numerics.Pi;
      use Ada.Numerics.Long_Elementary_Functions;
      Annual_Delay : constant := 13; -- months

      function Back_Step
        (Cur : Long_Float; RawV : Long_Float) return Long_Float
      is
         C : constant Long_Float := abs lagC;
         D : constant Long_Float := Cur - RawV;
         Pos : constant Long_Float := D + C;  -- candidate with Prev >= 0
         Neg : constant Long_Float := D - C;  -- candidate with Prev <  0
      begin
         -- If only one candidate is sign-consistent, take it.
         if Pos < 0.0 then
            return Neg;
         elsif Neg >= 0.0 then
            return Pos;
         end if;

         -- Ambiguous zone (both sign-consistent): pick a stable preference.
         -- Here: prefer keeping the same sign as the current value.
         if Cur > 0.0 then
            return Pos;
         elsif Cur < 0.0 then
            return Neg;
         else
            -- Cur == 0: pick the smaller magnitude to avoid a jump.
            if abs Pos <= abs Neg then
               return Pos;
            else
               return Neg;
            end if;
         end if;
      end Back_Step;

   begin
      for I in Raw'Range loop
         Start_Index := I;
         exit when Raw (I).Date > Start;
      end loop;
      --Text_IO.Put_Line(Raw(Start_Index).Date'Img & " " & Start_Index'Img);

      -- The first few values are sensitive to prior information so can adjust
      -- these for better fits in the early part of the time series.
      Res (Start_Index + 1).Value := 0.0;
      Res (Start_Index + 2).Value := iA; -- iB;
      Res (Start_Index + 3).Value := iA;
      for I in Start_Index + 4 .. Raw'Last - 1 loop
         if Modulation then
            Ramp :=
              Long_Float'Copy_Sign (lagC, Res (I - 1).Value) *
              (1.0 - Cos (2.0 * Pi * Raw (I).Date + iC));
         else
            Ramp := Long_Float'Copy_Sign (lagC, Res (I - 1).Value);
         end if;
         if I > Start_Index + Annual_Delay then
            FB2 := Raw (I - Annual_Delay).Value * iB;
         else
            FB2 := 0.0;
         end if;
         -- Res(I).Value := (Raw(I-1).Value + Raw(I).Value + Raw(I-2).Value + Raw(I+1).Value + Raw(I-3).Value)/5.0
         Res (I).Value :=
           (Raw (I).Value) / 1.0 +
           (lagA + mA * Cos (Pi * Raw (I).Date + mB)) * Res (I - 1).Value +
           lagB * Res (I - 2).Value - Ramp + FB2;
      end loop;

   -- If Start is inside the series, create "pre-history" by running backwards
   -- using inverted lag factors and inverted ramp sign.
      if Start_Index - 1 in Raw'Range and then Start_Index + 2 in Raw'Range
      then
         for I in reverse (Raw'First + 2) .. (Start_Index + 3) loop
            Res (I - 1).Value := Back_Step (Res (I).Value, Raw (I).Value);
         end loop;
      end if;

      return Res;
   end IIR_0;

   -- also called a boxcar filter, a moving average w/ window of 3 points
   function FIR
     (Raw : in Data_Pairs; Behind, Current, Ahead : in Long_Float)
      return Data_Pairs
   is
      Res : Data_Pairs := Raw;
   begin
      for I in Raw'First + 1 .. Raw'Last - 1 loop
         Res (I).Value :=
           Behind * Raw (I - 1).Value + Current * Raw (I).Value +
           Ahead * Raw (I + 1).Value;
      end loop;
      return Res;
   end FIR;

   -- Amplifies a tidal time series with an impulse array (i.e. Dirac comb)
   function Amplify
     (Raw : in Data_Pairs; Offset, Ramp, Start : in Long_Float)
      return Data_Pairs
   is
      Res : Data_Pairs := Raw;
   begin
      for I in Raw'Range loop
         -- Res(I).Value := Offset + Raw(I).Value * Impulse(Raw(I).Date + Ramp*(Raw(I).Date-Start));
         Res (I).Value :=
           Raw (I).Value *
           Impulse
             (Raw (I).Date + Ramp * (Raw (I).Date - Start) +
              Offset / 1_000_000.0 * (Raw (I).Date - Start)**2);
      end loop;
      return Res;
   end Amplify;

   -- Conventional tidal series summation or superposition of cycles
   function Tide_Sum_Diff
     (Template : in Data_Pairs; Constituents : in Long_Periods_Amp_Phase;
      Periods : in Long_Periods; Ref_Time : in Long_Float := 0.0;
      Scaling : in Long_Float := 1.0; Cos_Phase : in Boolean := True;
      Year_Len : in Long_Float := Year_Length; Integ : in Long_Float := 0.0;
      Ext_Forcing : in Data_Pairs := Empty_Data;
      Ext_Factor : in Long_Float := 0.0; Ext_Phase : in Long_Float := 0.0;
      Ext_Amp : in Long_Float := 0.0) return Data_Pairs
   is
      Pi : Long_Float := Ada.Numerics.Pi;
      Time : Long_Float;
      Res : Data_Pairs := Template;
      One : constant Long_Float := 1.0;
      Partition : constant Integer := 7; -- 4
      Pad : Integer;
      pragma Unreferenced (Integ);
   begin
      Pad := Integer (abs Ext_Factor * 1_000.0) mod Ealign;
      for I in Template'Range loop
         Time := Template (I).Date + Ref_Time;
         declare
            TF1, TF2 : Long_Float := 0.0;
            use Ada.Numerics.Long_Elementary_Functions;
         begin
            for J in Constituents'First .. Constituents'First + Partition loop
               declare
                  L : Amp_Phase renames Constituents (J);
                  Freq : Long_Float := Year_Len / Periods (J);
               begin
                  if Aliased_Period then
                     Freq := Long_Float'Remainder (Freq, One);
                  end if;
                  if Cos_Phase then
                     TF1 :=
                       TF1 +
                       L.Amplitude * (Cos (2.0 * Pi * Freq * Time + L.Phase));
                  else
                     TF1 :=
                       TF1 +
                       L.Amplitude * (Sin (2.0 * Pi * Freq * Time + L.Phase));
                  end if;
               end;
            end loop;
            for J in Constituents'First + Partition + 1 .. Constituents'Last
            loop
               declare
                  L : Amp_Phase renames Constituents (J);
                  Freq : Long_Float := Year_Len / Periods (J);
               begin
                  if Aliased_Period then
                     Freq := Long_Float'Remainder (Freq, One);
                  end if;
                  if Cos_Phase then
                     TF2 :=
                       TF2 +
                       L.Amplitude * (Cos (2.0 * Pi * Freq * Time + L.Phase));
                  else
                     TF2 :=
                       TF2 +
                       L.Amplitude * (Sin (2.0 * Pi * Freq * Time + L.Phase));
                  end if;
               end;
            end loop;
            Res (I) := (Time, TF1 + TF2 + Scaling * TF1 * TF2);
            if Ext_Forcing /= Empty_Data then
               for I in Res'First + Pad + 1 .. Res'Last - Pad loop
                  -- Res(I).Value := Res(I).Value + Ext_Factor * Ext_Forcing(I-2).Value + Ext_Phase * Ext_Forcing(I-1).Value;
                  --Res(I).Value := Res(I).Value + Ext_Factor * (Ext_Forcing(I-12).Value + Ext_Forcing(I-11).Value + Ext_Forcing(I-10).Value
                  --                                             + Ext_Forcing(I-9).Value + Ext_Forcing(I-8).Value + Ext_Forcing(I-7).Value
                  --                                             + Ext_Forcing(I-6).Value + Ext_Forcing(I-5).Value + Ext_Forcing(I-4).Value
                  --                                             + Ext_Forcing(I-3).Value + Ext_Forcing(I-2).Value + Ext_Forcing(I-1).Value)
                  --  + Ext_Phase * Ext_Forcing(I).Value;
                  if Ext_Factor > 0.0 then
                     Res (I).Value :=
                       Res (I).Value +
                       Ext_Phase * Ext_Forcing (I + Pad).Value +
                       Ext_Amp * Ext_Forcing (I + Pad - 1).Value;
                  else
                     Res (I).Value :=
                       Res (I).Value +
                       Ext_Phase * Ext_Forcing (I - Pad).Value +
                       Ext_Amp * Ext_Forcing (I - Pad - 1).Value;
                  end if;
               end loop;
            end if;
         end;
      end loop;
      return Res;
   end Tide_Sum_Diff;

-- Conventional tidal series summation or superposition of cycles
   function Tide_Sum_Custom
     (Template : in Data_Pairs; Constituents : in Long_Periods_Amp_Phase;
      Periods : in Long_Periods; Ref_Time : in Long_Float := 0.0;
      Scaling : in Long_Float := 1.0; Cos_Phase : in Boolean := True;
      Year_Len : in Long_Float := Year_Length; Integ : in Long_Float := 0.0;
      Ext_Forcing : in Data_Pairs := Empty_Data;
      Ext_Factor : in Long_Float := 0.0; Ext_Phase : in Long_Float := 0.0;
      Ext_Amp : in Long_Float := 0.0) return Data_Pairs
   is
      Pi : Long_Float := Ada.Numerics.Pi;
      Time : Long_Float;
      Res : Data_Pairs := Template;
      D, A, T, DD, TD, FF, DPlain, N, F : Long_Float;
      TF : Long_Float;
      pragma Unreferenced
        (Integ, Scaling, Cos_Phase, Ext_Forcing, Ext_Factor, Ext_Phase,
         Ext_Amp);
   begin
      for I in Template'Range loop
         Time := Template (I).Date + Ref_Time;
         declare -- uses 18.6y and 173
            use Ada.Numerics.Long_Elementary_Functions;
            L : Amp_Phase renames Constituents (2);
            Freq : Long_Float := Year_Len / Periods (2);
            HF : Amp_Phase renames Constituents (5);
            HP : Long_Float := Year_Len / Periods (5);
            NF : Amp_Phase renames Constituents (11);
            NP : Long_Float := Year_Len / Periods (11);
         begin
            D :=
              L.Amplitude *
              (Cos
                 (2.0 * Pi * Freq * Time + L.Phase +
                  HF.Amplitude * Cos (2.0 * Pi * HP * Time + HF.Phase) +
                  NF.Amplitude * Cos (2.0 * Pi * NP * Time + NF.Phase)));
         end;
         declare
            use Ada.Numerics.Long_Elementary_Functions;
            L : Amp_Phase renames Constituents (9);
            Freq : Long_Float := Year_Len / Periods (9);
            PF : Amp_Phase renames Constituents (14);
            PP : Long_Float := Year_Len / (2.0 * Periods (14));
            EF : Amp_Phase renames Constituents (15);  -- arbitrary
         begin --ampA*(ABS(SIN(freq*$A18+$E$13*ABS(SIN(freq*($A18+$E$14)))+phaseA))))
            A :=
              L.Amplitude *
              (Cos
                 (2.0 * Pi * Freq * Time + L.Phase +
                  PF.Amplitude *
                    abs
                    (Cos
                       (2.0 * Pi * PP * Time +
                        EF.Amplitude *
                          abs (Cos (2.0 * Pi * PP * Time + EF.Phase)) +
                        PF.Phase))));
         end;
         declare
            use Ada.Numerics.Long_Elementary_Functions;
            L : Amp_Phase renames Constituents (7);
            Freq : Long_Float := Year_Len / Periods (7);
         begin --ampA*(ABS(SIN(freq*$A18+$E$13*ABS(SIN(freq*($A18+$E$14)))+phaseA))))
            T := L.Amplitude * Cos (2.0 * Pi * Freq * Time + L.Phase);
         end;
         declare
            use Ada.Numerics.Long_Elementary_Functions;
            L : Amp_Phase renames Constituents (12);
            Freq : Long_Float := Year_Len / Periods (12);
         begin --ampA*(ABS(SIN(freq*$A18+$E$13*ABS(SIN(freq*($A18+$E$14)))+phaseA))))
            FF := L.Amplitude * Cos (2.0 * Pi * Freq * Time + L.Phase);
         end;
         declare
            use Ada.Numerics.Long_Elementary_Functions;
            L : Amp_Phase renames Constituents (6);
            Freq : Long_Float := Year_Len / Periods (6);
         begin --ampA*(ABS(SIN(freq*$A18+$E$13*ABS(SIN(freq*($A18+$E$14)))+phaseA))))
            TD := L.Amplitude * Cos (2.0 * Pi * Freq * Time + L.Phase);
         end;
         declare
            use Ada.Numerics.Long_Elementary_Functions;
            L : Amp_Phase renames Constituents (8);
            Freq : Long_Float := Year_Len / Periods (8);
         begin --ampA*(ABS(SIN(freq*$A18+$E$13*ABS(SIN(freq*($A18+$E$14)))+phaseA))))
            DD := L.Amplitude * Cos (2.0 * Pi * Freq * Time + L.Phase);
         end;
         declare
            use Ada.Numerics.Long_Elementary_Functions;
            L : Amp_Phase renames Constituents (1);
            Freq : Long_Float := Year_Len / Periods (2);
         begin
            DPlain := L.Amplitude * Cos (2.0 * Pi * Freq * Time + L.Phase);
         end;
         declare
            use Ada.Numerics.Long_Elementary_Functions;
            L : Amp_Phase renames Constituents (10);
            Freq : Long_Float := Year_Len / Periods (11);
         begin
            N := L.Amplitude * Cos (2.0 * Pi * Freq * Time + L.Phase);
         end;
         declare
            use Ada.Numerics.Long_Elementary_Functions;
            L : Amp_Phase renames Constituents (13);
            Freq : Long_Float := Year_Len / Periods (22);
         begin
            F := L.Amplitude * Cos (2.0 * Pi * Freq * Time + L.Phase);
         end;

         declare
            A1 : Long_Float renames Constituents (16).Amplitude;
            A2 : Long_Float renames Constituents (17).Amplitude;
            A3 : Long_Float renames Constituents (18).Amplitude;
            A4 : Long_Float renames Constituents (19).Amplitude;
            D1 : Long_Float renames Constituents (20).Amplitude;
            D2 : Long_Float renames Constituents (21).Amplitude;
            D3 : Long_Float renames Constituents (22).Amplitude;
            D4 : Long_Float renames Constituents (23).Amplitude;
            DA2 : Long_Float renames Constituents (24).Amplitude;
            D2A2 : Long_Float renames Constituents (25).Amplitude;
            D2A : Long_Float renames Constituents (26).Amplitude;
            DA3 : Long_Float renames Constituents (27).Amplitude;
            D3A : Long_Float renames Constituents (28).Amplitude;
            DA : Long_Float renames Constituents (29).Amplitude;
         begin
            TF :=
              F + N + DPlain + FF + TD + DD + T + A1 * A + A2 * A**2 +
              A3 * A**3 + A4 * A**4 + D1 * D + D2 * D**2 + D3 * D**3 +
              D4 * D**4 + DA2 * D * (A**2) + D2A2 * (D**2) * (A**2) +
              D2A * (D**2) * A + DA3 * D * (A**3) + D3A * (D**3) * A +
              DA * D * A;
            Res (I) := (Time, TF);
         end;
      end loop;
      return Res;
   end Tide_Sum_Custom;

   function Tide_Sum
     (Template : in Data_Pairs; Constituents : in Long_Periods_Amp_Phase;
      Periods : in Long_Periods; Ref_Time : in Long_Float := 0.0;
      Scaling : in Long_Float := 1.0; Cos_Phase : in Boolean := True;
      Year_Len : in Long_Float := Year_Length; Integ : in Long_Float := 0.0;
      Ext_Forcing : in Data_Pairs := Empty_Data;
      Ext_Factor : in Long_Float := 0.0; Ext_Phase : in Long_Float := 0.0;
      Ext_Amp : in Long_Float := 0.0) return Data_Pairs
   is
      Res : Data_Pairs := Template;
      pragma Unreferenced (Ref_Time);
   begin
      if Custom_Tide then
         Res :=
           Tide_Sum_Custom
             (Template, Constituents, Periods, 0.0, Scaling, Cos_Phase,
              Year_Len, Integ, Ext_Forcing, Ext_Factor, Ext_Phase, Ext_Amp);
      else
         Res :=
           Tide_Sum_Diff
             (Template, Constituents, Periods, 0.0, Scaling, Cos_Phase,
              Year_Len, Integ, Ext_Forcing, Ext_Factor, Ext_Phase, Ext_Amp);
      end if;
      return Res;
   end Tide_Sum;

   -- Stands for Laplace's Tidal Equation solution
   -- see Mathematical Geoenergy (2018), Ch.12
   function LTE
     (Forcing : in Data_Pairs; Wave_Numbers : in Modulations;
      Amp_Phase : in Modulations_Amp_Phase;
      Offset, K0, Trend, Accel : in Long_Float := 0.0;
      NonLin : in Long_Float := 1.0; Third : in Long_Float := 0.0)
      return Data_Pairs
   is
      Res : Data_Pairs := Forcing;
   begin
      for I in Forcing'Range loop
         declare
            LF : Long_Float := 0.0;
            use Ada.Numerics.Long_Elementary_Functions;
            Pi : Long_Float := Ada.Numerics.Pi;
         begin
            for J in Wave_Numbers'Range loop
               declare
                  M : GEM.LTE.Amp_Phase renames Amp_Phase (J);
                  SW : Long_Float;
               begin
                  SW :=
                    Sin
                      (2.0 * Pi * Wave_Numbers (J) * Res (I).Value + M.Phase) *
                    Exp (Res (I).Value * Third);
                  if SW < 0.0 then
                     SW := -(abs SW)**NonLin;
                  else
                     SW := SW**NonLin;
                  end if;
                  if Sinc > 0.0 then
                     LF :=
                       LF + M.Amplitude * SW / (abs (Res (I).Value) + Sinc);
                  else
                     LF := LF + M.Amplitude * SW;
                  end if;
               exception
                  when Constraint_Error =>
                     null;
               end;
            end loop;
            LF :=
              LF + Offset + K0 * Res (I).Value + Trend * Res (I).Date +
              Accel *
                (Res (I).Date -
                   Res (Forcing'First)
                     .Date)**2.0; -- K0 is wavenumber=0 solution
            Res (I).Value := LF;
         end;
      end loop;
      return Res;
   end LTE;

   -- some of the values may not be modifiesd so this is used to identify them
   function Is_Fixed (Value : in Long_Float) return Boolean is
      pragma Unreferenced (Value);
   begin
      return False;
   end Is_Fixed;

   --
   -- The following are utility functions for simulation
   --

   -- Pearson's Correlation Coefficient
   function CC (X, Y : in Data_Pairs) return Long_Float is
      N : Integer := 0; -- X'Length;
      sum_X, sum_Y, sum_XY, squareSum_X, squareSum_Y : Long_Float := 0.0;
      use Ada.Numerics.Long_Elementary_Functions;
      J : Integer; --:= Y'First;
      Denominator : Long_Float;
      Start, Stop : Integer;
   begin
      Start := X'First;
      Stop := X'Last;
      for I in X'Range loop
         if X (I).Value = 0.0 or Y (I).Value = 0.0 then
            Start := I + 1;
         else
            exit;
         end if;
      end loop;
      for I in reverse X'Range loop
         if X (I).Value = 0.0 or Y (I).Value = 0.0 then
            Stop := I - 1;
         else
            exit;
         end if;
      end loop;
      J := Y'First + (Start - X'First);
      for I in Start .. Stop loop
         -- sum of elements of array X.
         if X (I).Date > Start_Year and X (I).Date < End_Year and
           not (X (I).Value = 0.0 or Y (J).Value = 0.0)
         then

            sum_X := sum_X + X (I).Value;

            -- sum of elements of array Y.
            sum_Y := sum_Y + Y (J).Value;

            -- sum of X[i] * Y[i].
            sum_XY := sum_XY + X (I).Value * Y (J).Value;

            -- sum of square of array elements.
            squareSum_X := squareSum_X + X (I).Value * X (I).Value;
            squareSum_Y := squareSum_Y + Y (J).Value * Y (J).Value;
            N := N + 1;
         end if;

         J := J + 1;
      end loop;
      -- assert the denominator that it doesn't hit zero
      Denominator :=
        (Long_Float (N) * squareSum_X - sum_X * sum_X) *
        (Long_Float (N) * squareSum_Y - sum_Y * sum_Y);
      if Denominator <= 0.0 then
         return 0.0;
      else
         return (Long_Float (N) * sum_XY - sum_X * sum_Y) / Sqrt (Denominator);
      end if;
   exception
      when Constraint_Error =>
         return 0.0;
   end CC;

   -- A zero-crossing metric that is faster than CC which can be used to
   -- improve fits to time-series with many zero crossings, as the precise
   -- amplitude is not critical
   function Xing (X, Y : in Data_Pairs) return Long_Float is
      sum_absXY, sum_XY : Long_Float := 0.0;
      J : Integer := Y'First;
   begin
      for I in X'Range loop
         -- sum of X[i] * Y[i].
         sum_XY := sum_XY + X (I).Value * Y (J).Value;
         sum_absXY := sum_absXY + abs (X (J).Value) * abs (Y (J).Value);

         J := J + 1;
      end loop;
      return sum_XY / sum_absXY;
   end Xing;

   function Diff (Raw : in Data_Pairs) return Data_Pairs is
      Res : Data_Pairs := Raw;
   begin
      Res (Raw'First).Value := 0.0;
      for I in Raw'First + 1 .. Raw'Last loop
         Res (I).Value := Raw (I).Value - Raw (I - 1).Value;
      end loop;
      return Res;
   end Diff;

   function DER_CC (Model, Data : in Data_Pairs) return Long_Float is
   begin
      return CC (X => Diff (Model), Y => Diff (Data));
   end DER_CC;

   function RMS
     (X, Y : in Data_Pairs; Ref, Offset : in Long_Float) return Long_Float
   is
      sum_XY : Long_Float := 0.0;
      use Ada.Numerics.Long_Elementary_Functions;
      J : Integer := Y'First;
   begin
      for I in X'Range loop
         -- sum of (X[i] - Y[i])^2
         sum_XY := sum_XY + (X (I).Value - Y (J).Value + Offset)**2;
         J := J + 1;
      end loop;
      return 1.0 / (1.0 + Sqrt (sum_XY) / Ref);
   end RMS;

   function CID (X, Y : in Data_Pairs) return Long_Float is
      -- sum_XY,
      sum_X, sum_Y : Long_Float := 0.0;
      use Ada.Numerics.Long_Elementary_Functions;
   begin
      -- sum_XY := CC(X,Y);
      for I in X'First .. X'Last - 1 loop
         sum_X := sum_X + (X (I).Value - X (I + 1).Value)**2;
      end loop;
      for I in Y'First .. Y'Last - 1 loop
         sum_Y := sum_Y + (Y (I).Value - Y (I + 1).Value)**2;
      end loop;

      -- return sum_XY * Sqrt(Long_Float'Min(Sum_X, Sum_Y) / Long_Float'Max(Sum_X, Sum_Y));
      return
        Sqrt (Long_Float'Min (sum_X, sum_Y) / Long_Float'Max (sum_X, sum_Y));
   end CID;

   function DTW_Distance
     (X, Y : in Data_Pairs; Window_Size : Positive) return Long_Float
   is
      function Distance
        (X, Y : in Data_Pairs; Window_Size : Positive) return Long_Float
      is
         N : Positive := X'Length;
         type Real_Array is
           array (X'First - Window_Size .. X'Last + Window_Size) of Long_Float;
         DTW_Current, DTW_Previous : Real_Array := (others => Long_Float'Last);
      begin
         DTW_Previous (X'First) := 0.0;

         for I in X'First .. X'Last loop
            DTW_Current (X'First) := Long_Float'Last; -- Reset current row

            for J in
              Integer'Max (X'First, I - Window_Size) ..
                Integer'Min (X'Last, I + Window_Size)
            loop
               declare
                  Cost : Long_Float := abs (X (I).Value - Y (J).Value);
                  Min_Cost : Long_Float;
                  J_Previous : Integer := Integer'Max (J - 1, Y'First);
               begin
                  if Y (J).Value = 0.0 then
                     Cost := 0.0;
                  end if;
                  -- Compute minimum cost considering the DTW constraint
                  Min_Cost :=
                    Long_Float'Min
                      (Long_Float'Min
                         (DTW_Previous (J_Previous), DTW_Current (J_Previous)),
                       (if J > Y'First then DTW_Previous (J)
                        else Long_Float'Last));

                  DTW_Current (J) := Cost + Min_Cost;
               end;
               N := J;
            end loop;

            -- Swap the rows
            DTW_Previous := DTW_Current;
         end loop;

         return DTW_Current (N);
      end Distance;

      function Neg (X : in Data_Pairs) return Data_Pairs is
         N : Data_Pairs := X;
      begin
         for I in N'Range loop
            N (I).Value := -N (I).Value;
         end loop;
         return N;
      end Neg;
      Max : Long_Float;
   begin
      Max := Distance (Neg (Y), Y, Window_Size);
      return (Max - Distance (X, Y, Window_Size)) / Max;
   end DTW_Distance;

   procedure Sort (Arr : in out Data_Pairs) is
      Temp : Long_Float;
      Swapped : Boolean;
   begin
      for I in Arr'First .. Arr'Last - 1 loop
         Swapped := False;
         for J in Arr'First .. Arr'Last - I loop
            if Arr (J).Value > Arr (J + 1).Value then
               Temp := Arr (J).Value;
               Arr (J).Value := Arr (J + 1).Value;
               Arr (J + 1).Value := Temp;
               Swapped := True;
            end if;
         end loop;
         exit when not Swapped;
      end loop;
   end Sort;

   function EMD
     (X, Y : in Data_Pairs; Derivative : in Boolean := False) return Long_Float
   is

      function Wasserstein1 (Y1, Y2 : Data_Pairs) return Long_Float is
         Sorted_Y1 : Data_Pairs := Y1;
         Sorted_Y2 : Data_Pairs := Y2;
         Distance : Long_Float := 0.0;
      begin
         -- Sort both arrays
         Sort (Sorted_Y1);
         Sort (Sorted_Y2);

         -- Compute W1 distance (for equal-length arrays)
         for I in Y1'Range loop
            Distance :=
              Distance + abs (Sorted_Y1 (I).Value - Sorted_Y2 (I).Value);
         end loop;

         return Distance / Long_Float (Y1'Length);
      end Wasserstein1;

   begin
      if Derivative then
         return 1.0 / (1.0 + Wasserstein1 (Diff (X), Diff (Y)));
      else
         return 1.0 / (1.0 + Wasserstein1 (X, Y));
      end if;
   end EMD;

   function Median (Arr : in Data_Pairs) return Long_Float is
      Temp : Data_Pairs := Arr;
   begin
      Sort (Temp);
      return Temp ((Arr'Last - Arr'First) / 2 + Arr'First).Value;
   end Median;

   function Mean
     (Arr : in Data_Pairs; Absolute : in Boolean := False) return Long_Float
   is
      Temp : Long_Float := 0.0;
   begin
      for I in Arr'Range loop
         if Absolute then
            Temp := Temp + abs Arr (I).Value;
         else
            Temp := Temp + Arr (I).Value;
         end if;
      end loop;
      return Temp / Long_Float (Arr'Length);
   end Mean;

   function Scaled_Error_Metric (X, Y : in Data_Pairs) return Long_Float is

      -- Function to compute the Median Absolute Deviation (MAD)
      function Compute_MAD (Y : Data_Pairs) return Long_Float is
         Median_Y : constant Long_Float := Median (Y);
         Deviations : Data_Pairs := Y;
      begin
         for I in Y'Range loop
            Deviations (I).Value := abs (Y (I).Value - Median_Y);
         end loop;
         return Median (Deviations);
      end Compute_MAD;

      -- Function to compute the global-scale
      function Compute_Global_Scale
        (Y : Data_Pairs; Alpha : Long_Float := 0.01;
         Epsilon : Long_Float := 1.0e-10) return Long_Float
      is
         MAD_Y : constant Long_Float := Compute_MAD (Y);
         Mean_Abs_Y : constant Long_Float := Mean (Y, True);
         Scaled_Mean : constant Long_Float :=
           Alpha * (if Mean_Abs_Y > Epsilon then Mean_Abs_Y else Epsilon);
      begin
         return Long_Float'Max (MAD_Y, Scaled_Mean);
      end Compute_Global_Scale;

      -- Function to compute the Score
      function Compute_Score
        (X, Y : Data_Pairs; Alpha : Long_Float := 0.01) return Long_Float
      is
         Global_Scale : constant Long_Float := Compute_Global_Scale (Y, Alpha);
         Sum : Long_Float := 0.0;
      begin
         for I in Y'Range loop
            declare
               Error : constant Long_Float := X (I).Value - Y (I).Value;
               Scale : constant Long_Float :=
                 Long_Float'Max (Alpha * abs (Y (I).Value), Global_Scale);
            begin
               Sum := Sum + (Error**2) / (Scale**2);
            end;
         end loop;
         return Sum / Long_Float (Y'Length);
      end Compute_Score;

      Alpha : constant Long_Float := 0.01;
   begin
      return 1.0 / (1.0 + Compute_Score (X, Y, Alpha));
   end Scaled_Error_Metric;

   Pi : constant Long_Float := Ada.Numerics.Pi;
   Mult : constant Long_Float := GEM.Getenv ("FMULT", 1.006);
   Step : constant Long_Float := GEM.Getenv ("FSTEP", 0.18);
   F_Start : constant Long_Float := GEM.Getenv ("FSTART", 0.1);
   F_End : constant Long_Float := 1_000.0;

   function Min_Entropy_RMS (X, Y : in Data_Pairs) return Long_Float is
      use Ada.Numerics.Long_Elementary_Functions;
      Value, Sum : Long_Float := 0.0;
      S, C : Long_Float; -- cumulative
      F : Long_Float := F_Start;
      N : Integer := 1;
   begin
      loop
         S := 0.0;
         C := 0.0;
         for I in X'Range loop
            S := S + Sin (2.0 * Pi * F * X (I).Value) * Y (I).Value;
            C := C + Cos (2.0 * Pi * F * X (I).Value) * Y (I).Value;
         end loop;
         Value := Value + (S * S + C * C);
         Sum := Sum + Sqrt ((S * S + C * C));
         if Linear_Step then
            F := F + Step;
         else
            F := F * Mult;
         end if;
         exit when F > F_End;
         N := N + 1;
      end loop;
      Sum := Sum * Sum / Long_Float (N * N);
      Value := (Value - Sum) / Sum;
      --Text_IO.Put_Line(Value'Img & Sum'Img);
      return Value;
   end Min_Entropy_RMS;

   procedure ME_Power_Spectrum
     (Forcing, Model, Data : in Data_Pairs;
      Model_Spectrum, Data_Spectrum : out Data_Pairs; RMS : out Long_Float;
      Phase : in Boolean := False)
   is
      use Ada.Numerics.Long_Elementary_Functions;
      Model_S : Data_Pairs (Model'Range);
      Data_S : Data_Pairs (Data'Range);
      Value, Sum : Long_Float := 0.0;
      F : Long_Float := F_Start;
      S, C : Long_Float; -- cumulative
      N : Integer := 1;
   begin
      for J in Data'Range loop
         S := 0.0;
         C := 0.0;
         for I in Data'First + 8 .. Data'Last loop  -- remove Init value
            S := S + Sin (2.0 * Pi * F * Forcing (I).Value) * Data (I).Value;
            C := C + Cos (2.0 * Pi * F * Forcing (I).Value) * Data (I).Value;
         end loop;
         Data_S (J).Date := F;
         if Phase then
            Data_S (J).Value := S;
         else
            Data_S (J).Value := (S * S + C * C);
         end if;
         Value := Value + Data_S (J).Value;
         Sum := Sum + Sqrt (Data_S (J).Value);
         if Linear_Step then
            F := F + Step;
         else
            F := F * Mult;
         end if;
      end loop;
      F := F_Start;
      for J in Model'Range loop
         S := 0.0;
         C := 0.0;
         for I in Model'First + 8 .. Model'Last loop  -- remove Init value
            S := S + Sin (2.0 * Pi * F * Forcing (I).Value) * Model (I).Value;
            C := C + Cos (2.0 * Pi * F * Forcing (I).Value) * Model (I).Value;
         end loop;
         Model_S (J).Date := F;
         if Phase then
            Model_S (J).Value := S;
         else
            Model_S (J).Value := (S * S + C * C);
         end if;
         Value := Value + Model_S (J).Value;
         Sum := Sum + Sqrt (Model_S (J).Value);
         if Linear_Step then
            F := F + Step;
         else
            F := F * Mult;
         end if;
      end loop;
      Model_Spectrum := Model_S;
      Data_Spectrum := Data_S;
      Sum := Sum * Sum / Long_Float (N * N);
      RMS := (Value - Sum) / Sum;
   end ME_Power_Spectrum;

   function Min_Entropy_Power_Spectrum (X, Y : in Data_Pairs) return Long_Float
   is
      First : Positive := X'First;
      Last : Positive := X'Last;
      Mid : Positive := (First + Last) / 2;
      FD : Data_Pairs := Y (First .. Mid);
      LD : Data_Pairs := Y (Mid .. Last);
      RMS : Long_Float;
   begin
      if GEM.Getenv ("MERMS", False) then
         return Min_Entropy_RMS (X, Y);
      else
         ME_Power_Spectrum (X, FD, LD, FD, LD, RMS, False);
         return CC (Filter9Point (FD), Filter9Point (LD));
      end if;
   end Min_Entropy_Power_Spectrum;

   function FT_CC (Model, Data, Forcing : in Data_Pairs) return Long_Float is
      Model_S : Data_Pairs := Model;
      Data_S : Data_Pairs := Data;
      T : Data_Pairs := Forcing;
      RMS : Long_Float;
   begin
      for I in T'Range loop
         T (I).Value := T (I).Date;
      end loop;
      ME_Power_Spectrum
        (Forcing => T, Model => Model, Data => Data, Model_Spectrum => Model_S,
         Data_Spectrum => Data_S, RMS => RMS);
      Model_S := Window (Model_S, 2);
      Data_S := Window (Data_S, 2);
      -- return CC(Model_S, Data_S);
      return DTW_Distance (Model_S, Data_S, 9);
   end FT_CC;

   --
   procedure Dump (Model, Data : in Data_Pairs; Run_Time : Long_Float := 200.0)
   is
   begin
      for I in Data'Range loop
         Text_IO.Put_Line
           (Data (I).Date'Img & " " & Model (I).Value'Img & " " &
            Data (I).Value'Img);
         exit when Model (I).Date > Run_Time;
      end loop;
   end Dump;

   Running : Boolean := True;

   procedure Stop is
   begin
      Running := False;
   end Stop;

   function Halted return Boolean is
   begin
      return not Running;
   end Halted;

   procedure Continue is
   begin
      Running := True;
   end Continue;

   --
   -- protect the file from reentrancy
   --

   protected Safe is
      procedure Save
        (Model, Data, Forcing : in Data_Pairs; File_Name : in String);
      procedure Save
        (Model : in Data_Pairs; Mag : in Integer; File_Name : in String);
   end Safe;

   protected body Safe is
      procedure Save
        (Model, Data, Forcing : in Data_Pairs; File_Name : in String)
      is
         FT : Text_IO.File_Type;
         Model_S : Data_Pairs := Model;
         Data_S : Data_Pairs := Data;
         RMS : Long_Float;
      begin
         Text_IO.Create
           (File => FT, Name => File_Name, Mode => Text_IO.Out_File);
         if Model_S'Length < 10_000 then
            if Is_Minimum_Entropy then
               ME_Power_Spectrum
                 (Forcing => Model, Model => Forcing, Data => Data,
                  Model_Spectrum => Model_S, Data_Spectrum => Data_S,
                  RMS => RMS, Phase => False);
            else
               ME_Power_Spectrum
                 (Forcing => Forcing, Model => Model, Data => Data,
                  Model_Spectrum => Model_S, Data_Spectrum => Data_S,
                  RMS => RMS);
               Model_S := Window (Model_S, 2);
               Data_S := Window (Data_S, 2);
            end if;
         end if;
         for I in Data'Range loop
            Text_IO.Put_Line
              (FT,
               Data (I).Date'Img & ", " & Model (I).Value'Img & ", " &
               Data (I).Value'Img & ", " & Forcing (I).Value'Img & ", " &
               Data_S (I).Date'Img & ", " & Model_S (I).Value'Img & ", " &
               Data_S (I).Value'Img);
         end loop;
         Text_IO.Close (FT);
      end Save;

      procedure Save
        (Model : in Data_Pairs; Mag : in Integer; File_Name : in String)
      is
         FT : Text_IO.File_Type;
         Model_S : Data_Pairs := Expand (Model, Mag);
      begin
         if Mag > 1 then
            Text_IO.Create
              (File => FT, Name => File_Name, Mode => Text_IO.Out_File);
            for I in Model_S'Range loop
               Text_IO.Put_Line
                 (FT, Model_S (I).Date'Img & ", " & Model_S (I).Value'Img);
            end loop;
            Text_IO.Close (FT);
         end if;
      end Save;
   end Safe;

   procedure Save
     (Model, Data, Forcing : in Data_Pairs;
      File_Name : in String := "lte_results.csv")
   is
   begin
      Safe.Save (Model, Data, Forcing, File_Name);
      Safe.Save (Forcing, Magnify, "mag_" & File_Name);
   end Save;

   -- 3 point median
   function Median (Raw : in Data_Pairs) return Data_Pairs is
      Res : Data_Pairs := Raw;
   begin
      for I in Raw'First + 1 .. Raw'Last - 1 loop
         if Raw (I - 1).Value < Raw (I).Value then
            if Raw (I - 1).Value >= Raw (I + 1).Value then
               Res (I).Value := Raw (I - 1).Value;
            elsif Raw (I).Value < Raw (I + 1).Value then
               Res (I).Value := Raw (I).Value;
            else
               Res (I).Value := Raw (I + 1).Value;
            end if;
         else
            if Raw (I - 1).Value < Raw (I + 1).Value then
               Res (I).Value := Raw (I - 1).Value;
            elsif Raw (I).Value >= Raw (I + 1).Value then
               Res (I).Value := Raw (I).Value;
            else
               Res (I).Value := Raw (I + 1).Value;
            end if;
         end if;
      end loop;
      return Res;
   end Median;

   function Window
     (Raw : in Data_Pairs; Lobe_Width : in Positive) return Data_Pairs
   is
      Res : Data_Pairs := Raw;
   begin
      for I in Raw'First + Lobe_Width .. Raw'Last - Lobe_Width loop
         Res (I).Value := 0.0;
         for J in I - Lobe_Width .. I + Lobe_Width loop
            Res (I).Value := Res (I).Value + Raw (J).Value;
         end loop;
         Res (I).Value := Res (I).Value / Long_Float (2 * Lobe_Width + 1);
      end loop;
      return Res;
   end Window;

   -- Multiple Linear Regression types

   package MLR is new Ada.Numerics.Generic_Real_Arrays (Real => Long_Float);
   subtype Vector is MLR.Real_Vector;
   subtype Matrix is MLR.Real_Matrix;

   function To_Matrix
     (Source : Vector; Column_Vector : Boolean := True) return Matrix
   is
      Result : Matrix (1 .. 1, Source'Range);
   begin
      for Column in Source'Range loop
         Result (1, Column) := Source (Column);
      end loop;
      if Column_Vector then
         return MLR.Transpose (Result);
      else
         return Result;
      end if;
   end To_Matrix;

   function To_Row_Vector
     (Source : Matrix; Column : Positive := 1) return Vector
   is
      Result : Vector (Source'Range (1));
   begin
      for Row in Result'Range loop
         Result (Row) := Source (Row, Column);
      end loop;
      return Result;
   end To_Row_Vector;

   function Regression_Coefficients
     (Source : Vector; Regressors : Matrix) return Vector
   is
      Result : Matrix (Regressors'Range (2), 1 .. 1);
      Nil : Vector (1 .. 0);
   begin
      if Source'Length /= Regressors'Length (1) then
         raise Constraint_Error;
      end if;
      declare
         Regressors_T : constant Matrix := MLR.Transpose (Regressors);
         use MLR;
      begin
         Result :=
           MLR.Inverse (Regressors_T * Regressors) * Regressors_T *
           To_Matrix (Source);
      end;
      return To_Row_Vector (Source => Result);
   exception
      when Constraint_Error => -- Singular
         Text_IO.Put_Line ("Singular result, converging?");
         delay 1.0;
         return Nil; -- Source;  -- Doesn't matter, will give a junk result
   end Regression_Coefficients;

   procedure Regression_Factors
     (Data_Records : in Data_Pairs;
      --First, Last,  -- Training Interval

      NM : in Positive; -- # modulations
      Forcing : in Data_Pairs;  -- Value @ Time
      -- Factors_Matrix : in out Matrix;

      DBLT : in Periods; DALTAP : out Amp_Phases; DALEVEL : out Long_Float;
      DAK0 : out Long_Float; Secular_Trend : in out Long_Float;
      Accel : out Long_Float; Singular : out Boolean;
      Third : in Long_Float := 0.0)
   is

      use Ada.Numerics.Long_Elementary_Functions;
      First : Integer := Data_Records'First;
      Last : Integer := Data_Records'Last;
      Pi : Long_Float := Ada.Numerics.Pi;
      Trend : Boolean := Secular_Trend > 0.0;
      -- Add_Trend : Integer := Integer(Secular_Trend);
      Add_Trend : Integer := 2 * Boolean'Pos (Trend);
      Num_Coefficients : constant Integer :=
        2 + NM * 2 + Add_Trend; -- !!! sin + cos mod
      RData : Vector (1 .. Last - First + 1);
      Factors_Matrix : Matrix (1 .. Last - First + 1, 1 .. Num_Coefficients);
      Value : Long_Float;
   begin
      for I in First .. Last loop
         RData (I - First + 1) := Data_Records (I).Value;
         Factors_Matrix (I - First + 1, 1) := 1.0;  -- DC offset
         Factors_Matrix (I - First + 1, 2) := Forcing (I).Value;
         for K in DBLT'First .. NM loop  -- D.B.LT'First = 1
            Value :=
              Sin (2.0 * Pi * DBLT (K) * Forcing (I).Value) *
              Exp (Forcing (I).Value * Third);
            Factors_Matrix (I - First + 1, 3 + (K - 1) * 2) := Value;
            Value :=
              Cos (2.0 * Pi * DBLT (K) * Forcing (I).Value) *
              Exp (Forcing (I).Value * Third);
            Factors_Matrix (I - First + 1, 4 + (K - 1) * 2) := Value;
         end loop;
         if Trend then
            Factors_Matrix (I - First + 1, Num_Coefficients - 1) :=
              Forcing (I).Date;
            Factors_Matrix (I - First + 1, Num_Coefficients) :=
              (Forcing (I).Date - Forcing (First).Date)**2.0;
         end if;
      end loop;

      declare
         Coefficients : constant Vector := -- MLR.
           Regression_Coefficients
             (Source => RData, Regressors => Factors_Matrix);
         K : Integer := 4;
      begin
         if Coefficients'Length = 0 then
            Singular := True;
         else
            DALEVEL := Coefficients (1);
            DAK0 := Coefficients (2);
            for I in 1 .. NM loop  -- if odd
               DALTAP (K / 2 - 1).Amplitude :=
                 Sqrt
                   (Coefficients (K - 1) * Coefficients (K - 1) +
                    Coefficients (K) * Coefficients (K));
               DALTAP (K / 2 - 1).Phase :=
                 Arctan (Coefficients (K), Coefficients (K - 1));
               K := K + 2;
            end loop;
            if Trend then
               Secular_Trend := Coefficients (Num_Coefficients - 1); --!!!
               Accel := Coefficients (Num_Coefficients); --!!!
            else
               Secular_Trend := 0.0;
               Accel := 0.0;
            end if;
            Singular := False;
         end if;
      end;
   end Regression_Factors;

   procedure Put
     (Value : in Long_Float; Text : in String := "";
      New_Line : in Boolean := False)
   is
   begin
      Ada.Long_Float_Text_IO.Put (Value, Fore => 4, Aft => 11, Exp => 0);
      Text_IO.Put (Text);
      if New_Line then
         Text_IO.New_Line;
      end if;
   end Put;

   function Filter9Point (Raw : in Data_Pairs) return Data_Pairs is
      Res : Data_Pairs := Raw;
   begin
      for I in Raw'First + 1 .. Raw'Last - 1 loop -- 4
         Res (I).Value :=
           0.25 * (+Raw (I - 1).Value) + 0.5 * (+Raw (I).Value) +
           0.25 * (+Raw (I + 1).Value);
      end loop;
      return Res;
   end Filter9Point;

end GEM.LTE.Primitives;
