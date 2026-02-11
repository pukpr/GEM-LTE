--  ============================================================================
--  GEM - General Environment Management Utilities
--  ============================================================================
--
--  PURPOSE:
--    Provides centralized configuration management for the GEM-LTE system.
--    Implements a three-level configuration hierarchy:
--      1. Response files (.resp) - Key-value pairs for batch configuration
--      2. Environment variables - Runtime overrides
--      3. Command-line flags - Boolean switches
--
--  KEY FEATURES:
--    - Unified Getenv() interface with type overloading (String, Integer,
--      Long_Float, Boolean) - single API for all configuration access
--    - Response file parsing: Reads <executable>.resp for parameter lists
--    - String parsing utilities: Convert strings to integer/float arrays
--    - Command-line option detection
--
--  CONFIGURATION PRECEDENCE (highest to lowest):
--    1. Command-line flags (e.g., -p, -w, -l) - treated as TRUE booleans
--    2. Response file parameters - loaded at startup
--    3. Environment variables - system-level configuration
--    4. Default values - provided in Getenv() calls
--
--  RESPONSE FILE FORMAT:
--    <executable>.resp contains name-value pairs:
--      CLIMATE_INDEX nino4.dat
--      TIMEOUT 0.0
--      TRIGGER 0.0
--      MAXLOOPS 100000
--      ... etc
--
--  USAGE EXAMPLE:
--    Timeout : Long_Float := GEM.Getenv("TIMEOUT", 0.999_999);
--    MaxLoops : Integer := GEM.Getenv("MAXLOOPS", 100_000);
--    UseTrend : Boolean := GEM.Getenv("TREND", False);
--
--  ============================================================================

with Ada.Environment_Variables;
with Ada.Command_Line.Response_File;
with Text_IO;
with Ada.Integer_Text_IO;
with Ada.Long_Float_Text_IO;
with Ada.Directories;

package body GEM is

   --  Enumeration of all recognized configuration parameters
   --  Used for response file parsing validation
   type Options is
     (ALIAS, EXCLUDE, ALTERNATE, CLIMATE_INDEX, DLOD_REF, EVERY, FILTER, FLIP,
      FORCING, FSTEP, IMPA, IMPB, IMPC, IMPD, IMPULSE, IR, MAXH, MAXLOOPS,
      MERMS, METRIC, MLR, NH, NM, PARETO, RESET, SAMPLING, SCALING, SINPOW,
      SPLIT_LOW, SPLIT_TRAINING, SPREAD_CYCLE, SPREAD_MIN, SPREAD_MAX,
      STARTING_METRIC, STEP, THRESHOLD, THRESHOLD_ACTION, TRAIN_END,
      TRAIN_START, TREND, YEAR, F9, DECAY, FMULT, FSTART, LOCKF, R2M, SYM,
      YTRIM, NONLIN, IDATE, TEST_ONLY, TIMEOUT, NUMBER_OF_PROCESSORS);

   type Option_Pair is record
      Name, Value : Ada.Command_Line.Response_File.String_Access;
   end record;

   type Options_list is array (Options) of Option_Pair;
   OL : Options_list;  -- Parsed response file parameters

   --  TODO: Can remove - Original idea was to directly rename Environment_Variables.Value
   --  but this prevented debugging and didn't support response file/command-line precedence.
   --  Custom implementation below provides 3-level hierarchy (cmd-line > resp > env).
   --function Getenv (Name : in String; Default : in String) return String renames
   --  Ada.Environment_Variables.Value; -- works but can't debug

   --  =========================================================================
   --  Getenv: Unified configuration retrieval with 3-level hierarchy
   --
   --  PRECEDENCE:
   --    1. Command-line flags (boolean only): "-p" returns "TRUE"
   --    2. Response file parameters: Parsed from <executable>.resp
   --    3. Environment variables: System-level configuration
   --    4. Default value: Provided by caller
   --
   --  Returns String representation, caller converts via 'Value attribute
   --  =========================================================================
   function Getenv (Name : in String; Default : in String) return String is
      Str : constant String := Ada.Environment_Variables.Value (Name, Default);
      use type Ada.Command_Line.Response_File.String_Access;
   begin
      -- Look for command line options that may be booleans
      for I in 1 .. Ada.Command_Line.Argument_Count loop
         if Ada.Command_Line.Argument (I) = Name then
            return "TRUE";
         end if;
      end loop;
      -- Look for any parameters parsed from response file
      for EV in Options loop
         if OL (EV).Name /= null and then OL (EV).Name.all = Name then
            return OL (EV).Value.all;
         end if;
      end loop;
      -- Otherwise return the Env Var
      return Str;
   end Getenv;

   --  Type-specific overloads: Convert string result to appropriate type
   function Getenv (Name : in String; Default : in Integer) return Integer is
   begin
      return Integer'Value (Getenv (Name, Integer'Image (Default)));
   end Getenv;

   function Getenv (Name : in String; Default : in Float) return Float is
   begin
      return Float'Value (Getenv (Name, Float'Image (Default)));
   end Getenv;

   function Getenv
     (Name : in String; Default : in Long_Integer) return Long_Integer
   is
   begin
      return Long_Integer'Value (Getenv (Name, Long_Integer'Image (Default)));
   end Getenv;

   function Getenv
     (Name : in String; Default : in Long_Float) return Long_Float
   is
   begin
      return Long_Float'Value (Getenv (Name, Long_Float'Image (Default)));
   end Getenv;

   function Getenv (Name : in String; Default : in Boolean) return Boolean is
   begin
      return Boolean'Value (Getenv (Name, Boolean'Image (Default)));
   end Getenv;

   --  Set environment variable at runtime
   procedure Setenv (Name : in String; Value : in String) is
   begin
      Ada.Environment_Variables.Set (Name, Value);
   end Setenv;

   --  Clear environment variable
   procedure Clear (Name : in String) is
   begin
      Ada.Environment_Variables.Clear (Name);
   end Clear;

   --  =========================================================================
   --  String Parsing Utilities
   --  =========================================================================

   --  Parse space-separated integers from string into array
   --  Used for NH (harmonics) parameter: "2 3 5 7" -> [2,3,5,7]
   function S_to_I (S : in String) return Ns is
      use Ada.Integer_Text_IO;
      List : Ns (1 .. 100); -- magic number
      N, L : Integer := 0;
      Index : Integer := 1;
   begin
      loop
         Get (S (L + 1 .. S'Last), N, L);
         List (Index) := N;
         Index := Index + 1;
      end loop;
   exception
      when others =>
         return List (1 .. Index - 1);
   end S_to_I;

   --  Parse space-separated floats from string into array
   --  Similar to S_to_I but for Long_Float values
   function S_to_LF (S : in String) return Fs is
      use Ada.Long_Float_Text_IO;
      List : Fs (1 .. 100); -- magic number
      N : Long_Float;
      L : Integer := 0;
      Index : Integer := 1;
   begin
      loop
         Get (S (L + 1 .. S'Last), N, L);
         List (Index) := N;
         Index := Index + 1;
      end loop;
   exception
      when others =>
         return List (1 .. Index - 1);
   end S_to_LF;

   --  Check if command-line option flag exists
   --  Used for -p (print), -w (write), -l (load legacy) flags
   function Command_Line_Option_Exists (Option : in String) return Boolean is
   begin
      for I in 1 .. Ada.Command_Line.Argument_Count loop
         if Ada.Command_Line.Argument (I) = Option then
            Text_IO.Put_Line ("EXITING");
            return True;
         end if;
      end loop;
      return False;
   end Command_Line_Option_Exists;

   --  =========================================================================
   --  Read_Response_File: Parse configuration from <executable>.resp
   --
   --  ALGORITHM:
   --    1. Construct filename from executable name + ".resp"
   --    2. Parse into argument list (name-value pairs)
   --    3. Match names against Options enumeration
   --    4. Store in OL array for Getenv() lookup
   --    5. Print configuration summary to console
   --
   --  FORMAT:
   --    Each line: OPTION_NAME value
   --    Example: CLIMATE_INDEX nino4.dat
   --  =========================================================================
   procedure Read_Response_File is
      FN : constant String :=
        Ada.Directories.Simple_Name (Ada.Command_Line.Command_Name) & ".resp";
      use type Ada.Command_Line.Response_File.String_Access;
   begin
      Text_IO.Put_Line ("RESPONSE FILE");
      declare
         L : Ada.Command_Line.Response_File.Argument_List :=
           Ada.Command_Line.Response_File.Arguments_From (FN);
      begin
         --  TODO: Can remove - Debug output for response file parsing was useful
         --  during development but now redundant (lines 207-209 print all values).
         --Ada.Text_IO.Put_Line(L(2*I-1).all & " " & L(2*I).all);
         for I in L'First .. L'Last / 2 loop
            for EV in Options loop
               if EV'Image = L (2 * I - 1).all then
                  OL (EV).Name := L (2 * I - 1);
                  OL (EV).Value := L (2 * I);
               end if;
            end loop;
         end loop;
      end;
      for EV in Options loop
         if OL (EV).Name = null then
            Text_IO.Put_Line (EV'Image & " not set, use default");
         else
            Text_IO.Put_Line (OL (EV).Name.all & " " & OL (EV).Value.all);
         end if;
      end loop;

   exception
      when Ada.Command_Line.Response_File.File_Does_Not_Exist =>
         Text_IO.Put_Line (FN & " not found");
   end Read_Response_File;

begin

   --  Package initialization: Load response file at startup
   Read_Response_File;

end GEM;
