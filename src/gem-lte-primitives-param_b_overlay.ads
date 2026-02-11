--  ============================================================================
--  GEM.LTE.Primitives.Param_B_Overlay - Safe Array Overlay for Param_B
--  ============================================================================
--
--  PURPOSE:
--    Provides type-safe wrapper around memory overlay used for random descent
--    optimization. Encapsulates the unsafe 'Address clause with runtime checks.
--
--  BACKGROUND:
--    The random descent algorithm (Walker.Markov) needs to randomly modify
--    parameters without knowledge of the Param_B record structure. Rather than
--    using unsafe C-style memory overlay, this package provides:
--    - Named constants for field positions in overlay array
--    - Runtime verification of layout assumptions
--    - Type-safe interface for creating overlays
--
--  FIELD INDICES IN OVERLAY ARRAY:
--    Scalars (18 fields):
--      1. Offset    7. DelB      13. Year    
--      2. bg        8. Asym      14. mA      
--      3. ImpA      9. Ann1      15. mP      
--      4. ImpB     10. Ann2      16. shiftT  
--      5. ImpC     11. Sem1      17. init    
--      6. DelA     12. Sem2      18. IR      
--
--    Arrays (variable size):
--      19..19+NLP*2-1: LPAP (each constituent = 2 floats: Amplitude, Phase)
--      19+NLP*2..end:  LT modulations
--
--  USAGE:
--    Size := Param_B_Overlay.Overlay_Size (D.NLP, D.NLT);
--    Set : LF_Array (1 .. Size);
--    for Set'Address use D.B.Offset'Address;
--    Param_B_Overlay.Verify_Layout (D.B, Set);  -- Runtime check
--
--  ============================================================================

with GEM.LTE.Primitives.Shared;

package GEM.LTE.Primitives.Param_B_Overlay is

   use GEM.LTE.Primitives.Shared;
   
   --  Named constants for scalar field positions in overlay array
   Offset_Index  : constant := 1;
   BG_Index      : constant := 2;
   ImpA_Index    : constant := 3;
   ImpB_Index    : constant := 4;
   ImpC_Index    : constant := 5;
   DelA_Index    : constant := 6;
   DelB_Index    : constant := 7;
   Asym_Index    : constant := 8;
   Ann1_Index    : constant := 9;
   Ann2_Index    : constant := 10;
   Sem1_Index    : constant := 11;
   Sem2_Index    : constant := 12;
   IR_Index      : constant := 13;
   Year_Index    : constant := 14;
   MA_Index      : constant := 15;
   MP_Index      : constant := 16;
   ShiftT_Index  : constant := 17;
   Init_Index    : constant := 18;
   
   Scalar_Field_Count : constant := 18;
   
   --  Calculate overlay array size for given discriminants
   --  Formula: Scalar_Fields + (NLP * 2) + NLT
   function Overlay_Size (NLP, NLT : Integer) return Positive;
   
   --  Verify that overlay array correctly maps to Param_B record
   --  Raises Constraint_Error if layout assumptions are violated
   --  This should be called once at initialization for safety
   procedure Verify_Layout (
      P : in Param_B;
      Expected_Size : in Positive);
   
   --  Get index of first LPAP element (Amplitude of constituent 1)
   --  LPAP indices: First_LPAP_Index .. First_LPAP_Index + NLP*2 - 1
   function First_LPAP_Index return Positive;
   
   --  Get index of first LT element
   --  LT indices: First_LT_Index(NLP) .. First_LT_Index(NLP) + NLT - 1
   function First_LT_Index (NLP : Integer) return Positive;
   
   --  Get index for specific LPAP constituent's amplitude
   --  Constituent 1..NLP, returns index for Amplitude field
   function LPAP_Amplitude_Index (Constituent : Positive) return Positive;
   
   --  Get index for specific LPAP constituent's phase
   --  Constituent 1..NLP, returns index for Phase field
   function LPAP_Phase_Index (Constituent : Positive) return Positive;

end GEM.LTE.Primitives.Param_B_Overlay;
