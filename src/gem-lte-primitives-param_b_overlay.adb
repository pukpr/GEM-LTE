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

end GEM.LTE.Primitives.Param_B_Overlay;
