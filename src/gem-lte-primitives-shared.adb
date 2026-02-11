--  ============================================================================
--  GEM.LTE.Primitives.Shared - Parameter Storage and I/O Management
--  ============================================================================
--
--  PURPOSE:
--    Manages the shared parameter state for multi-threaded GEM-LTE optimization.
--    Provides thread-safe storage, persistence (read/write), and serialization
--    (text .par files and JSON) of model parameters.
--
--  KEY RESPONSIBILITIES:
--    1. Protected Server: Thread-safe parameter sharing across worker threads
--    2. Parameter Persistence: Save/load from .par files (text format)
--    3. JSON Serialization: Load parameters from JSON format
--    4. Parameter Validation: Read with format checking and error handling
--
--  DATA STRUCTURES:
--    - Param_S: Main parameter record (defined in .ads)
--      - B: Variable parameters (offsets, impulses, LTE coefficients)
--      - A: Constant parameters (tidal periods, amplitudes, phases)
--      - C: Harmonic configuration array
--    - Protected Server: Manages single shared instance with mutex access
--
--  FILE FORMATS:
--    - .par files: Human-readable text format (4-char tags + values)
--      Example: "offs  0.12345", "impA  1.23456", "ltep  2.71828"
--    - .json files: JSON format for easier configuration management
--      Example: {"offs": 0.12345, "impA": 1.23456, "ltep": [2.71828, ...]}
--    - .parms files: Binary Direct_IO format (legacy, fast but not portable)
--
--  NAMING CONVENTION:
--    Files are named after the executable:
--    - enso_opt.exe.par (general parameters)
--    - enso_opt.exe.nino4.dat.par (climate-index-specific parameters)
--    - enso_opt.exe.json (JSON format)
--
--  THREAD SAFETY:
--    Server protected object ensures only one thread can modify shared
--    parameters at a time. Get() blocks until Put() has been called at
--    least once.
--
--  ============================================================================

with Ada.Direct_IO;
with Ada.Command_Line;
with Ada.Text_IO;
with GNAT.OS_Lib;
with Ada.Long_Float_Text_IO;
with Ada.Integer_Text_IO;
with GNATCOLL.JSON;
with Ada.Strings.Unbounded;
with Ada.Characters.Latin_1;
with Ada.IO_Exceptions;
with Ada.Directories;

package body GEM.LTE.Primitives.Shared is

   CI : constant String := GEM.Getenv (Name => "CLIMATE_INDEX", Default => "");

   --  Access type required for protected object storage
   --  (cannot directly store discriminated record in protected object)
   type Param_P is access all Param_S;

   --  =========================================================================
   --  Server: Protected object for thread-safe parameter sharing
   --
   --  PATTERN:
   --    Single-writer-multiple-reader with blocking Get(). First Put()
   --    allocates storage and sets Available flag, subsequent Put() calls
   --    reuse the allocated space. Get() blocks until Available is true.
   --
   --  USAGE:
   --    - Put() called by optimization thread when better parameters found
   --    - Get() called by worker threads to initialize or update their state
   --    - One-time allocation on first Put() avoids repeated heap allocation
   --  =========================================================================
   protected Server is
      procedure Put (P : in Param_S);
      entry Get (P : in out Param_S);
   private
      Available : Boolean := False;
      Params : Param_P;
   end Server;

   ------------
   -- Server --
   ------------

   protected body Server is

      --  Store parameters (allocates on first call, updates thereafter)
      procedure Put (P : in Param_S) is
      begin
         if not Available then -- only do once to allocate
            Params := new Param_S'(P);
         else  -- use allocated space
            Params.all := P;
         end if;
         Available := True;
      end Put;

      --  Retrieve parameters (blocks until first Put() completes)
      entry Get (P : in out Param_S) when Available is
      begin
         P := Params.all;
      end Get;

   end Server;

   --  Public wrappers for Server access
   procedure Put (P : in Param_S) is
   begin
      Server.Put (P);
   end Put;

   function Get (N_Tides, N_Modulations : in Integer) return Param_S is
      P : Param_S (N_Tides, N_Modulations);
   begin
      Server.Get (P);
      return P;
   end Get;

   --  =========================================================================
   --  Output Formatting Helpers
   --  =========================================================================

   --  Write parameter name + Long_Float value in fixed format
   procedure Put
     (FT : in Ada.Text_IO.File_Type; Text : in String; Value : in Long_Float)
   is
   begin
      Ada.Text_IO.Put (FT, Text & " ");
      Ada.Long_Float_Text_IO.Put (FT, Value, Fore => 4, Aft => 11, Exp => 0);
      Ada.Text_IO.New_Line (FT);
   end Put;

   --  Write parameter name + Integer value
   procedure Put
     (FT : in Ada.Text_IO.File_Type; Text : in String; Value : in Integer)
   is
   begin
      Ada.Text_IO.Put (FT, Text & " " & Integer'Image (Value));
      Ada.Text_IO.New_Line (FT);
   end Put;

   --  Write triplet (period, amplitude, phase) for tidal constituents
   procedure Put (FT : Ada.Text_IO.File_Type; V1, V2, V3 : in Long_Float) is
   begin
      Ada.Long_Float_Text_IO.Put (FT, V1, Fore => 4, Aft => 11, Exp => 0);
      Ada.Long_Float_Text_IO.Put (FT, V2, Fore => 4, Aft => 11, Exp => 0);
      Ada.Long_Float_Text_IO.Put (FT, V3, Fore => 4, Aft => 11, Exp => 0);
      Ada.Text_IO.New_Line (FT);
   end Put;

   --  =========================================================================
   --  Write: Save parameters to text .par files
   --
   --  Generates two files:
   --    1. <executable>.par - General parameters for all climate indices
   --    2. <executable>.<climate_index>.par - Index-specific parameters
   --
   --  Format: 4-character tag followed by space and value
   --  Example: "offs  0.12345678901"
   --  =========================================================================
   procedure Write (D : in Param_S) is
      FN : constant String :=
        Ada.Directories.Simple_Name (Ada.Command_Line.Command_Name) & ".par";
      FN2 : constant String :=
        Ada.Directories.Simple_Name (Ada.Command_Line.Command_Name) & "." &
        CI & ".par";
      FT : Ada.Text_IO.File_Type;
   begin
      Ada.Text_IO.Create (FT, Ada.Text_IO.Out_File, FN);
      Put (FT, "offs", D.B.Offset);
      Put (FT, "bg  ", D.B.bg);
      Put (FT, "impA", D.B.ImpA);
      Put (FT, "impB", D.B.ImpB);
      Put (FT, "impC", D.B.ImpC);
      Put (FT, "delA", D.B.DelA);
      Put (FT, "delB", D.B.DelB);
      Put (FT, "asym", D.B.Asym);
      Put (FT, "ann1", D.B.Ann1);
      Put (FT, "ann2", D.B.Ann2);
      Put (FT, "sem1", D.B.Sem1);
      Put (FT, "sem2", D.B.Sem2);
      Put (FT, "year", D.B.Year);
      Put (FT, "IR  ", D.B.IR);
      Put (FT, "ma  ", D.B.mA);
      Put (FT, "mp  ", D.B.mP);
      Put (FT, "shfT", D.B.shiftT);
      Put (FT, "init", D.B.init);
      for I in D.B.LPAP'Range loop
         Ada.Text_IO.Put (FT, " ");
         Put (FT, D.A.LP (I), D.B.LPAP (I).Amplitude, D.B.LPAP (I).Phase);
      end loop;
      for I in D.B.LT'Range loop
         Put (FT, "ltep", D.B.LT (I));
      end loop;
      for I in D.C'Range loop
         Put (FT, "harm", D.C (I));
         exit when D.C (I) = 0;
      end loop;
      Ada.Text_IO.Close (FT);
      -- per data set
      Ada.Text_IO.Create (FT, Ada.Text_IO.Out_File, FN2);
      Put (FT, "offs", D.B.Offset);
      Put (FT, "bg  ", D.B.bg);
      Put (FT, "impA", D.B.ImpA);
      Put (FT, "impB", D.B.ImpB);
      Put (FT, "impC", D.B.ImpC);
      Put (FT, "delA", D.B.DelA);
      Put (FT, "delB", D.B.DelB);
      Put (FT, "asym", D.B.Asym);
      Put (FT, "ann1", D.B.Ann1);
      Put (FT, "ann2", D.B.Ann2);
      Put (FT, "sem1", D.B.Sem1);
      Put (FT, "sem2", D.B.Sem2);
      Put (FT, "year", D.B.Year);
      Put (FT, "IR  ", D.B.IR);
      Put (FT, "ma  ", D.B.mA);
      Put (FT, "mp  ", D.B.mP);
      Put (FT, "shfT", D.B.shiftT);
      Put (FT, "init", D.B.init);
      for I in D.B.LT'Range loop
         Put (FT, "ltep", D.B.LT (I));
      end loop;
      for I in D.C'Range loop
         Put (FT, "harm", D.C (I));
         exit when D.C (I) = 0;
      end loop;
      Ada.Text_IO.Close (FT);
   end Write;

   --  =========================================================================
   --  Save: Persist parameters to binary .parms file and text .par files
   --
   --  Uses Ada.Direct_IO for fast binary serialization (legacy format).
   --  Also writes human-readable .par files via Write() procedure.
   --  =========================================================================
   procedure Save (P : in Param_S) is
      subtype PS is Param_S (P.NLP, P.NLT);
      package DIO is new Ada.Direct_IO (PS);
      FN : constant String :=
        Ada.Directories.Simple_Name (Ada.Command_Line.Command_Name) & ".parms";
      use DIO;
      FT : File_Type;
   begin
      Create (FT, Out_File, FN);
      Write (FT, P);
      Close (FT);
      --  TODO: Can remove - Command-line flag 'r' was intended for conditional
      --  writing but Write() is now always called. The flag check is redundant.
      --if GEM.Command_Line_Option_Exists("r") then
      Write (P);
      --end if;
      -- Server.Put (P);
   end Save;

   --  =========================================================================
   --  Input Parsing Helpers - Read parameters from text files
   --  =========================================================================

   --  Read 4-char tag + Long_Float value, validate tag matches expected
   procedure Read
     (FT : in Ada.Text_IO.File_Type; Text : in String;
      Value : in out Long_Float)
   is
      Str : String (1 .. 100);
      N : Integer;
   begin
      Ada.Text_IO.Get_Line (FT, Str, N);
      if Text = Str (1 .. 4) then
         Ada.Long_Float_Text_IO.Get (Str (5 .. N), Value, N);
      else
         Ada.Text_IO.Put_Line (Text & " mismatches " & Str (1 .. N));
         GNAT.OS_Lib.OS_Exit (0);
      end if;
   end Read;

   --  Read 4-char tag + Integer value, validate tag matches expected
   procedure Read
     (FT : in Ada.Text_IO.File_Type; Text : in String; Value : in out Integer)
   is
      Str : String (1 .. 100);
      N : Integer;
   begin
      Ada.Text_IO.Get_Line (FT, Str, N);
      if Text = Str (1 .. 4) then
         Ada.Integer_Text_IO.Get (Str (5 .. N), Value, N);
      else
         Ada.Text_IO.Put_Line (Text & " mismatches " & Str (1 .. N));
         GNAT.OS_Lib.OS_Exit (0);
      end if;
   end Read;

   procedure Read (FT : in Ada.Text_IO.File_Type; V1, V2, V3 : out Long_Float)
   is
      Str : String (1 .. 100);
      N : Integer;
      Floats : Fs (1 .. 3);
   begin
      Ada.Text_IO.Get_Line (FT, Str, N);
      Floats := S_to_LF (Str (1 .. N));
      V1 := Floats (1);
      V2 := Floats (2);
      V3 := Floats (3);
   exception
      when others =>
         Ada.Text_IO.Put_Line ("Finished " & Str (1 .. N));
   end Read;

   function File_To_String (Name : in String) return String is
      FT : Ada.Text_IO.File_Type;
      Line : String (1 .. 1_024);
      Last : Natural;
      Contents : Ada.Strings.Unbounded.Unbounded_String :=
        Ada.Strings.Unbounded.Null_Unbounded_String;
   begin
      Ada.Text_IO.Open (FT, Ada.Text_IO.In_File, Name);
      begin
         while not Ada.Text_IO.End_Of_File (FT) loop
            Ada.Text_IO.Get_Line (FT, Line, Last);
            Ada.Strings.Unbounded.Append (Contents, Line (1 .. Last));
            Ada.Strings.Unbounded.Append (Contents, Ada.Characters.Latin_1.LF);
         end loop;
         Ada.Text_IO.Close (FT);
      exception
         when others =>
            if Ada.Text_IO.Is_Open (FT) then
               Ada.Text_IO.Close (FT);
            end if;
            raise;
      end;
      return Ada.Strings.Unbounded.To_String (Contents);
   end File_To_String;

   function Get_Number (Value : in GNATCOLL.JSON.JSON_Value) return Long_Float
   is
      use GNATCOLL.JSON;
   begin
      case Kind (Value) is
         when JSON_Float_Type =>
            return Get_Long_Float (Value);
         when JSON_Int_Type =>
            declare
               Int_Value : constant Long_Long_Integer := Get (Value);
            begin
               return Long_Float (Int_Value);
            end;
         when others =>
            Ada.Text_IO.Put_Line ("JSON number error:" & Kind (Value)'Img);
            GNAT.OS_Lib.OS_Exit (0);
      end case;
   end Get_Number;

   function Get_Integer (Value : in GNATCOLL.JSON.JSON_Value) return Integer is
      use GNATCOLL.JSON;
   begin
      case Kind (Value) is
         when JSON_Int_Type =>
            declare
               Int_Value : constant Long_Long_Integer := Get (Value);
            begin
               return Integer (Int_Value);
            end;
         when JSON_Float_Type =>
            return Integer (Get_Long_Float (Value));
         when others =>
            Ada.Text_IO.Put_Line ("JSON integer error:" & Kind (Value)'Img);
            GNAT.OS_Lib.OS_Exit (0);
      end case;
   end Get_Integer;

   function Get_Float_Field
     (Data : in GNATCOLL.JSON.JSON_Value; Name : in String;
      Default : in Long_Float) return Long_Float
   is
      use GNATCOLL.JSON;
   begin
      if Kind (Data) = JSON_Object_Type and then Has_Field (Data, Name) then
         return Get_Number (Get (Data, Name));
      end if;
      return Default;
   end Get_Float_Field;

   procedure Read_JSON_Float_Array
     (Data : in GNATCOLL.JSON.JSON_Value; Name : in String;
      Target : in out Modulations)
   is
      use GNATCOLL.JSON;
      Arr : JSON_Array;
      Count : Natural;
   begin
      if Kind (Data) /= JSON_Object_Type or else not Has_Field (Data, Name)
      then
         return;
      end if;
      Arr := Get (Data, Name);
      Count := Length (Arr);
      declare
         J : Positive := 1;
      begin
         for I in Target'Range loop
            exit when J > Count;
            Target (I) := Get_Number (Get (Arr, J));
            J := J + 1;
         end loop;
      end;
   exception
      when others =>
         Ada.Text_IO.Put_Line ("JSON array error " & Name);
         GNAT.OS_Lib.OS_Exit (0);
   end Read_JSON_Float_Array;

   procedure Read_JSON_Int_Array
     (Data : in GNATCOLL.JSON.JSON_Value; Name : in String; Target : in out Ns)
   is
      use GNATCOLL.JSON;
      Arr : JSON_Array;
      Count : Natural;
   begin
      if Kind (Data) /= JSON_Object_Type or else not Has_Field (Data, Name)
      then
         return;
      end if;
      Arr := Get (Data, Name);
      Count := Length (Arr);
      declare
         J : Positive := 1;
      begin
         for I in Target'Range loop
            exit when J > Count;
            Target (I) := Get_Integer (Get (Arr, J));
            J := J + 1;
         end loop;
      end;
   exception
      when others =>
         Ada.Text_IO.Put_Line ("JSON array error " & Name);
         GNAT.OS_Lib.OS_Exit (0);
   end Read_JSON_Int_Array;

   --  =========================================================================
   --  JSON Reading Functions - Modern configuration format
   --
   --  Provides flexible JSON-based parameter loading with fallback to
   --  text .par format. JSON format is easier to edit and validate.
   --  =========================================================================

   --  Read JSON array of Long_Float triplets (period, amplitude, phase)
   --  Used for LPAP (Long Period Amplitude Phase) tidal constituent data
   procedure Read_JSON_LPAP
     (Data : in GNATCOLL.JSON.JSON_Value; Name : in String; D : in out Param_S)
   is
      use GNATCOLL.JSON;
      Arr : JSON_Array;
      Count : Natural;
      procedure Apply (Period, Amp, Phase : in Long_Float) is
      begin
         for I in D.B.LPAP'Range loop
            if abs (GEM.LTE.LP (I)) > 0.0
              and then abs (Period - GEM.LTE.LP (I)) <=
                abs (GEM.LTE.LP (I)) * 0.01
            then
               D.B.LPAP (I).Amplitude := Amp;
               D.B.LPAP (I).Phase := Phase;
               exit;
            end if;
         end loop;
      end Apply;
   begin
      if Kind (Data) /= JSON_Object_Type or else not Has_Field (Data, Name)
      then
         return;
      end if;
      Arr := Get (Data, Name);
      Count := Length (Arr);
      if Count = 0 then
         return;
      end if;
      if Kind (Get (Arr, 1)) = JSON_Array_Type then
         for I in 1 .. Count loop
            declare
               Triplet : constant JSON_Array := Get (Get (Arr, I));
            begin
               if Length (Triplet) >= 3 then
                  Apply
                    (Get_Number (Get (Triplet, 1)),
                     Get_Number (Get (Triplet, 2)),
                     Get_Number (Get (Triplet, 3)));
               end if;
            end;
         end loop;
      else
         declare
            Index : Positive := 1;
         begin
            while Index + 2 <= Count loop
               Apply
                 (Get_Number (Get (Arr, Index)),
                  Get_Number (Get (Arr, Index + 1)),
                  Get_Number (Get (Arr, Index + 2)));
               Index := Index + 3;
            end loop;
         end;
      end if;
   exception
      when others =>
         Ada.Text_IO.Put_Line ("JSON array error " & Name);
         GNAT.OS_Lib.OS_Exit (0);
   end Read_JSON_LPAP;

   --  Read complete parameter set from JSON file
   --  Returns True on success, False if file not found or parse error
   --  Include_LPAP flag controls whether tidal constituent data is loaded
   function Read_JSON
     (Name : in String; D : in out Param_S; Include_LPAP : in Boolean)
      return Boolean
   is
      use GNATCOLL.JSON;
      Result : Read_Result;
      Data : JSON_Value;
   begin
      declare
         Text : constant String := File_To_String (Name);
      begin
         Result := Read (Text);
      exception
         when Ada.Text_IO.Name_Error | Ada.IO_Exceptions.Name_Error =>
            return False;
         when others =>
            return False;
      end;
      if not Result.Success then
         Ada.Text_IO.Put_Line (Format_Parsing_Error (Result.Error));
         GNAT.OS_Lib.OS_Exit (0);
      end if;
      Data := Result.Value;
      D.B.Offset := Get_Float_Field (Data, "offs", D.B.Offset);
      D.B.bg := Get_Float_Field (Data, "bg", D.B.bg);
      D.B.ImpA := Get_Float_Field (Data, "impA", D.B.ImpA);
      D.B.ImpB := Get_Float_Field (Data, "impB", D.B.ImpB);
      D.B.ImpC := Get_Float_Field (Data, "impC", D.B.ImpC);
      D.B.DelA := Get_Float_Field (Data, "delA", D.B.DelA);
      D.B.DelB := Get_Float_Field (Data, "delB", D.B.DelB);
      D.B.Asym := Get_Float_Field (Data, "asym", D.B.Asym);
      D.B.Ann1 := Get_Float_Field (Data, "ann1", D.B.Ann1);
      D.B.Ann2 := Get_Float_Field (Data, "ann2", D.B.Ann2);
      D.B.Sem1 := Get_Float_Field (Data, "sem1", D.B.Sem1);
      D.B.Sem2 := Get_Float_Field (Data, "sem2", D.B.Sem2);
      D.B.Year := Get_Float_Field (Data, "year", D.B.Year);
      D.B.IR :=
        Get_Float_Field (Data, "IR", Get_Float_Field (Data, "ir", D.B.IR));
      D.B.mA := Get_Float_Field (Data, "ma", D.B.mA);
      D.B.mP := Get_Float_Field (Data, "mp", D.B.mP);
      D.B.shiftT :=
        Get_Float_Field
          (Data, "shfT", Get_Float_Field (Data, "shft", D.B.shiftT));
      D.B.init := Get_Float_Field (Data, "init", D.B.init);
      if Include_LPAP then
         Read_JSON_LPAP (Data, "LPAP", D);
         Read_JSON_LPAP (Data, "lpap", D);
      end if;
      Read_JSON_Float_Array (Data, "ltep", D.B.LT);
      Read_JSON_Float_Array (Data, "LT", D.B.LT);
      Read_JSON_Float_Array (Data, "lt", D.B.LT);
      Read_JSON_Int_Array (Data, "harm", D.C);
      Read_JSON_Int_Array (Data, "harmonics", D.C);
      for I in D.A.LP'Range loop
         D.A.LP (I) := GEM.LTE.LP (I);
      end loop;
      return True;
   exception
      when others =>
         Ada.Text_IO.Put_Line ("JSON read error " & Name);
         return False;
         -- GNAT.OS_Lib.Os_Exit(0);
   end Read_JSON;

   --  =========================================================================
   --  Read: Load parameters from text .par files (fallback to JSON)
   --
   --  Tries JSON first, falls back to text .par format. Loads:
   --    1. General parameters from <executable>.par or .json
   --    2. Index-specific overrides from <executable>.<climate_index>.par
   --
   --  Validates tag names and exits on mismatch to catch file corruption.
   --  =========================================================================
   procedure Read (D : in out Param_S) is
      Exec : constant String :=
        Ada.Directories.Simple_Name (Ada.Command_Line.Command_Name);
      Base : constant String :=
        (if Exec'Length > 0 and then Exec (Exec'Last) = '.' then
           Exec (Exec'First .. Exec'Last - 1)
         else Exec);
      FN : constant String := Base & ".par";
      FN2 : constant String := Base & "." & CI & ".par";
      FN_JSON : constant String := Base & ".json";
      FN2_JSON : constant String := Base & "." & CI & ".json";
      FT : Ada.Text_IO.File_Type;
   begin
      if not Read_JSON (FN_JSON, D, True) then
         Ada.Text_IO.Open (FT, Ada.Text_IO.In_File, FN);
         Read (FT, "offs", D.B.Offset);
         Read (FT, "bg  ", D.B.bg);
         Read (FT, "impA", D.B.ImpA);
         Read (FT, "impB", D.B.ImpB);
         Read (FT, "impC", D.B.ImpC);
         Read (FT, "delA", D.B.DelA);
         Read (FT, "delB", D.B.DelB);
         Read (FT, "asym", D.B.Asym);
         Read (FT, "ann1", D.B.Ann1);
         Read (FT, "ann2", D.B.Ann2);
         Read (FT, "sem1", D.B.Sem1);
         Read (FT, "sem2", D.B.Sem2);
         Read (FT, "year", D.B.Year);
         Read (FT, "IR  ", D.B.IR);
         Read (FT, "ma  ", D.B.mA);
         Read (FT, "mp  ", D.B.mP);
         Read (FT, "shfT", D.B.shiftT);
         Read (FT, "init", D.B.init);
         for I in D.B.LPAP'Range loop
            Read (FT, D.A.LP (I), D.B.LPAP (I).Amplitude, D.B.LPAP (I).Phase);
            D.A.LP (I) := GEM.LTE.LP (I); -- OVERRIDE!
         end loop;
         for I in D.B.LT'Range loop
            Read (FT, "ltep", D.B.LT (I));
         end loop;

         begin
            for I in D.C'Range loop
               Read (FT, "harm", D.C (I));
            end loop;
            Ada.Text_IO.Close (FT);
         exception
            when Ada.Text_IO.End_Error =>
               Ada.Text_IO.Put_Line ("Closing " & FN);
               if Ada.Text_IO.Is_Open (FT) then
                  Ada.Text_IO.Close (FT);
               end if;
            when others =>
               if Ada.Text_IO.Is_Open (FT) then
                  Ada.Text_IO.Close (FT);
               end if;
               raise;
         end;
      end if;

      if not Read_JSON (FN2_JSON, D, False) then
         Ada.Text_IO.Open (FT, Ada.Text_IO.In_File, FN2);
         Read (FT, "offs", D.B.Offset);
         Read (FT, "bg  ", D.B.bg);
         Read (FT, "impA", D.B.ImpA);
         Read (FT, "impB", D.B.ImpB);
         Read (FT, "impC", D.B.ImpC);
         Read (FT, "delA", D.B.DelA);
         Read (FT, "delB", D.B.DelB);
         Read (FT, "asym", D.B.Asym);
         Read (FT, "ann1", D.B.Ann1);
         Read (FT, "ann2", D.B.Ann2);
         Read (FT, "sem1", D.B.Sem1);
         Read (FT, "sem2", D.B.Sem2);
         Read (FT, "year", D.B.Year);
         Read (FT, "IR  ", D.B.IR);
         Read (FT, "ma  ", D.B.mA);
         Read (FT, "mp  ", D.B.mP);
         Read (FT, "shfT", D.B.shiftT);
         Read (FT, "init", D.B.init);
         for I in D.B.LT'Range loop
            Read (FT, "ltep", D.B.LT (I));
         end loop;
         begin
            for I in D.C'Range loop
               Read (FT, "harm", D.C (I));
            end loop;
            Ada.Text_IO.Close (FT);
         exception
            when Ada.Text_IO.End_Error =>
               Ada.Text_IO.Put_Line ("Closing " & FN2);
               if Ada.Text_IO.Is_Open (FT) then
                  Ada.Text_IO.Close (FT);
               end if;
            when others =>
               if Ada.Text_IO.Is_Open (FT) then
                  Ada.Text_IO.Close (FT);
               end if;
               Ada.Text_IO.Put_Line ("? Opening " & FN2);
         end;
      end if;
   exception
      when others =>
         if Ada.Text_IO.Is_Open (FT) then
            Ada.Text_IO.Close (FT);
         end if;
         Ada.Text_IO.Put_Line ("? Opening " & FN2);
   end Read;

   --  =========================================================================
   --  Load: Initialize parameters from file or legacy binary format
   --
   --  Command-line flags control behavior:
   --    -l : Load from legacy .parms binary file
   --    -p : Dump parameters and exit
   --    -w : Write parameters to .par files and exit
   --  =========================================================================
   procedure Load (P : in out Param_S) is
      subtype PS is Param_S (P.NLP, P.NLT);
      package DIO is new Ada.Direct_IO (PS);
      FN : constant String :=
        Ada.Directories.Simple_Name (Ada.Command_Line.Command_Name) & ".parms";
      use DIO;
      FT : File_Type;
   begin
      if GEM.Command_Line_Option_Exists ("l") then  -- legacy file
         begin
            Open (FT, In_File, FN);
            Read (FT, P);
            Close (FT);
         exception
            when Name_Error =>
               Ada.Text_IO.Put_Line ("No PARMS file: " & FN);
            when others =>
               Ada.Text_IO.Put_Line ("Error:" & FN);
               if Is_Open (FT) then
                  Close (FT);
               end if;
         end;
      else
         Read (P);
      end if;

      if GEM.Command_Line_Option_Exists ("p") then
         Dump (P);
         GNAT.OS_Lib.OS_Exit (0);
      elsif GEM.Command_Line_Option_Exists ("w") then
         Write (P);
         GNAT.OS_Lib.OS_Exit (0);
         --  TODO: Can remove - Command-line flag 'r' for Read was planned but
         --  never implemented. Read() is called automatically in normal flow.
         --elsif GEM.Command_Line_Option_Exists("r") then
         --   Read(P);
      end if;
   end Load;

   --  =========================================================================
   --  Dump: Print all parameters in human-readable markdown format
   --
   --  Outputs parameters with percentage differences from reference values.
   --  Used with -p command-line flag for parameter inspection.
   --  =========================================================================
   procedure Dump (D : in Param_S) is
      function Percent (A, B : in Long_Float) return String is
      begin
         return ", " & Integer'Image (Integer ((B - A) / A * 100.0));
      end Percent;
   begin
      Ada.Text_IO.Put_Line ("```");
      Put (D.B.Offset, " :offset:", NL);
      Put (D.B.bg, " :bg:", NL);
      Put (D.B.ImpA, " :impA:", NL);
      Put (D.B.ImpB, " :impB:", NL);
      Put (D.B.ImpC, " :impC:", NL);
      Put (D.B.DelA, ":delA:", NL);
      Put (D.B.DelB, ":delB:", NL);
      Put (D.B.Asym, ":asym:", NL);
      Put (D.B.Ann1, ":ann1:", NL);
      Put (D.B.Ann2, ":ann2:", NL);
      Put (D.B.Sem1, ":sem1:", NL);
      Put (D.B.Sem2, ":sem2:", NL);
      Put (D.B.Year, ":year:", NL);
      Put (D.B.IR, ":IR:", NL);
      Put (D.B.mA, " :mA:", NL);
      Put (D.B.mP, " :mP:", NL);
      Put (D.B.shiftT, " :shiftT:", NL);
      Put (D.B.init, " :init:", NL);
      Ada.Text_IO.Put_Line ("---- Tidal ----");
      for I in D.B.LPAP'Range loop
         Put (D.A.LP (I), ", ");
         Put (D.B.LPAP (I).Amplitude, ", ");
         Put
           (D.B.LPAP (I).Phase,
            ", " & I'Img &
            Percent (GEM.LTE.LPRef (I).Amplitude, D.B.LPAP (I).Amplitude) &
            ", " & GEM.LTE.LPRef (I).Amplitude'Img,
            NL);
      end loop;
      --  TODO: Can remove - Debug output for LTE static parameters was used
      --  during development to verify parameter loading. Now redundant as main
      --  Dump output shows all critical parameters.
      --Ada.Text_IO.Put_Line("---- LTE static ----");
      --for I in D.B.LT'Range loop
      --   Put(D.B.LT(I), "", NL);
      --end loop;
   end Dump;

end GEM.LTE.Primitives.Shared;
