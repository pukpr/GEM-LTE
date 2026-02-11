--  ============================================================================
--  GEM.LTE.Primitives.Param_B_Overlay - Implementation
--  ============================================================================

with Text_IO;

package body GEM.LTE.Primitives.Param_B_Overlay is

   function Overlay_Size (NLP, NLT : Integer) return Positive is
   begin
      return Scalar_Field_Count + (NLP * 2) + NLT;
   end Overlay_Size;
   
   function First_LPAP_Index return Positive is
   begin
      return Scalar_Field_Count + 1;  -- 19
   end First_LPAP_Index;
   
   function First_LT_Index (NLP : Integer) return Positive is
   begin
      return Scalar_Field_Count + (NLP * 2) + 1;
   end First_LT_Index;
   
   function LPAP_Amplitude_Index (Constituent : Positive) return Positive is
   begin
      --  Each constituent has 2 floats: Amplitude, Phase
      --  Constituent 1: indices 19 (Amp), 20 (Phase)
      --  Constituent 2: indices 21 (Amp), 22 (Phase)
      --  Formula: First_LPAP_Index + (Constituent - 1) * 2
      return First_LPAP_Index + (Constituent - 1) * 2;
   end LPAP_Amplitude_Index;
   
   function LPAP_Phase_Index (Constituent : Positive) return Positive is
   begin
      --  Phase comes after Amplitude
      return LPAP_Amplitude_Index (Constituent) + 1;
   end LPAP_Phase_Index;
   
   procedure Verify_Layout (
      P : in Param_B;
      Expected_Size : in Positive)
   is
      --  Calculate size from record's bit size
      --  Subtract discriminant storage (implementation-dependent)
      Calculated_Size : constant Positive := P'Size / Long_Float'Size - 1;
      
      Manual_Size : constant Positive := Overlay_Size (P.NLP, P.NLT);
   begin
      Text_IO.Put_Line ("=== Param_B Overlay Verification ===");
      Text_IO.Put_Line ("  NLP (tidal constituents) = " & P.NLP'Image);
      Text_IO.Put_Line ("  NLT (modulations) = " & P.NLT'Image);
      Text_IO.Put_Line ("  Expected overlay size = " & Expected_Size'Image);
      Text_IO.Put_Line ("  Calculated (P'Size / LF'Size - 1) = " & Calculated_Size'Image);
      Text_IO.Put_Line ("  Manual (18 + NLP*2 + NLT) = " & Manual_Size'Image);
      
      if Calculated_Size /= Expected_Size then
         Text_IO.Put_Line ("  ERROR: Calculated size mismatch!");
         raise Constraint_Error with 
            "Param_B overlay size mismatch: expected " & Expected_Size'Image &
            ", calculated " & Calculated_Size'Image;
      end if;
      
      if Manual_Size /= Expected_Size then
         Text_IO.Put_Line ("  ERROR: Manual size mismatch!");
         raise Constraint_Error with
            "Param_B overlay size mismatch: expected " & Expected_Size'Image &
            ", manual " & Manual_Size'Image;
      end if;
      
      --  Verify field indices are in range
      if LPAP_Phase_Index (P.NLP) > Expected_Size then
         raise Constraint_Error with
            "LPAP array extends beyond overlay bounds";
      end if;
      
      if First_LT_Index (P.NLP) + P.NLT - 1 > Expected_Size then
         raise Constraint_Error with
            "LT array extends beyond overlay bounds";
      end if;
      
      Text_IO.Put_Line ("  ✓ Layout verification PASSED");
      Text_IO.Put_Line ("  Scalar fields: 1.." & Scalar_Field_Count'Image);
      Text_IO.Put_Line ("  LPAP array: " & First_LPAP_Index'Image & 
                       ".." & Integer'Image(First_LPAP_Index + P.NLP * 2 - 1));
      Text_IO.Put_Line ("  LT array: " & Integer'Image(First_LT_Index(P.NLP)) & 
                       ".." & Integer'Image(First_LT_Index(P.NLP) + P.NLT - 1));
      Text_IO.Put_Line ("");
   end Verify_Layout;
   
   procedure Debug_Print_Field_Mapping (
      P : in Param_B;
      Set : in LF_Array)
   is
      use Text_IO;
      
      --  Check if debug mode is enabled
      Debug_Mode : constant String := GEM.Getenv ("OVERLAY_DEBUG", "0");
   begin
      if Debug_Mode /= "1" then
         return;  -- Debug mode not enabled
      end if;
      
      Put_Line ("=== OVERLAY DEBUG: Field-by-Field Mapping ===");
      Put_Line ("");
      Put_Line ("Scalar Fields (Direct Record Access vs Overlay Array):");
      Put_Line ("  Index  Field      Record Value      Overlay Value     Match?");
      Put_Line ("  -----  ---------  ----------------  ----------------  ------");
      
      --  Verify each scalar field
      declare
         type Field_Check is record
            Index : Positive;
            Name  : String (1 .. 9);
            Value : Long_Float;
         end record;
         
         Fields : constant array (1 .. Scalar_Field_Count) of Field_Check := (
            (Offset_Index, "Offset   ", P.Offset),
            (BG_Index,     "bg       ", P.bg),
            (ImpA_Index,   "ImpA     ", P.ImpA),
            (ImpB_Index,   "ImpB     ", P.ImpB),
            (ImpC_Index,   "ImpC     ", P.ImpC),
            (DelA_Index,   "DelA     ", P.DelA),
            (DelB_Index,   "DelB     ", P.DelB),
            (Asym_Index,   "Asym     ", P.Asym),
            (Ann1_Index,   "Ann1     ", P.Ann1),
            (Ann2_Index,   "Ann2     ", P.Ann2),
            (Sem1_Index,   "Sem1     ", P.Sem1),
            (Sem2_Index,   "Sem2     ", P.Sem2),
            (IR_Index,     "IR       ", P.IR),
            (Year_Index,   "Year     ", P.Year),
            (MA_Index,     "mA       ", P.mA),
            (MP_Index,     "mP       ", P.mP),
            (ShiftT_Index, "shiftT   ", P.shiftT),
            (Init_Index,   "init     ", P.init)
         );
         
         Matches : Natural := 0;
      begin
         for F of Fields loop
            declare
               Record_Val : constant Long_Float := F.Value;
               Overlay_Val : constant Long_Float := Set (F.Index);
               Match : constant Boolean := (Record_Val = Overlay_Val);
            begin
               Put ("  " & F.Index'Image);
               Put ("    " & F.Name);
               Put ("  " & Record_Val'Image);
               Put ("  " & Overlay_Val'Image);
               Put_Line ("  " & (if Match then "✓" else "✗"));
               
               if Match then
                  Matches := Matches + 1;
               end if;
            end;
         end loop;
         
         Put_Line ("");
         Put_Line ("Scalar fields matched: " & Matches'Image & 
                   " / " & Scalar_Field_Count'Image);
      end;
      
      --  Show LPAP array mapping
      Put_Line ("");
      Put_Line ("LPAP Array (Tidal Constituents):");
      Put_Line ("  Const  Amplitude_Idx  Phase_Idx  Amp_Value     Phase_Value");
      Put_Line ("  -----  -------------  ---------  ------------  ------------");
      
      for I in 1 .. P.NLP loop
         declare
            Amp_Idx : constant Positive := LPAP_Amplitude_Index (I);
            Pha_Idx : constant Positive := LPAP_Phase_Index (I);
         begin
            Put ("  " & I'Image);
            Put ("    " & Amp_Idx'Image);
            Put ("           " & Pha_Idx'Image);
            Put ("       " & Set(Amp_Idx)'Image);
            Put_Line ("  " & Set(Pha_Idx)'Image);
         end;
      end loop;
      
      --  Show LT array mapping
      Put_Line ("");
      Put_Line ("LT Array (Modulation Periods):");
      Put_Line ("  Index  Value");
      Put_Line ("  -----  ----------------");
      
      declare
        LT_Start : constant Positive := First_LT_Index (P.NLP);
      begin
         for I in 0 .. P.NLT - 1 loop
            declare
               Idx : constant Positive := LT_Start + I;
            begin
               Put_Line ("  " & Idx'Image & "   " & Set(Idx)'Image);
            end;
         end loop;
      end;
      
      Put_Line ("");
      Put_Line ("=== End of Overlay Debug ===");
      Put_Line ("");
   end Debug_Print_Field_Mapping;

end GEM.LTE.Primitives.Param_B_Overlay;
