with Ada.Environment_Variables;
with Ada.Command_Line.Response_File;
with Text_IO;
with Ada.Integer_Text_IO;
with Ada.Long_Float_Text_IO;
with Ada.Directories;

package body GEM is

   type Options is (ALIAS, EXCLUDE, ALTERNATE, CLIMATE_INDEX, DLOD_REF, EVERY, FILTER, FLIP, FORCING,
                    FSTEP, IMPA, IMPB, IMPC, IMPD, IMPULSE, IR, MAXH, MAXLOOPS,
                    MERMS, METRIC, MLR, NH, NM, PARETO, RESET, SAMPLING,
                    SCALING, SINPOW, SPLIT_LOW, SPLIT_TRAINING, SPREAD_CYCLE,
                    SPREAD_MIN, SPREAD_MAX, STARTING_METRIC, STEP, THRESHOLD,
                    THRESHOLD_ACTION, TRAIN_END, TRAIN_START, TREND, YEAR, F9,
                    DECAY, FMULT, FSTART, LOCKF, R2M, SYM, YTRIM, NONLIN, IDATE, TEST_ONLY, TIMEOUT, NUMBER_OF_PROCESSORS
                   );

   type Option_Pair is
      record
         Name, Value : Ada.Command_Line.Response_File.String_Access;
      end record;

   type Options_list is array(Options) of Option_Pair;
   OL : Options_List;

   --function Getenv (Name : in String; Default : in String) return String renames
   --  Ada.Environment_Variables.Value; -- works but can't debug

   function Getenv (Name : in String; Default : in String) return String is
      Str : constant String := Ada.Environment_Variables.Value(Name,Default);
      use type Ada.Command_Line.Response_File.String_Access;
   begin
      -- Look for command line options that may be booleans
      for I in 1 .. Ada.Command_Line.Argument_Count loop
         if Ada.Command_Line.Argument(I) = Name then
            return "TRUE";
         end if;
      end loop;
      -- Look for any parameters parsed from response file
      for EV in Options loop
         if OL(EV).Name /= null and then
           OL(EV).Name.all = Name then
            return OL(EV).Value.all;
         end if;
      end loop;
      -- Otherwise return the Env Var
      return Str;
   end Getenv;

   function Getenv (Name : in String; Default : in Integer) return Integer is
   begin
      return Integer'Value (Getenv (Name, Integer'Image (Default)));
   end Getenv;

   function Getenv (Name : in String; Default : in Float) return Float is
   begin
      return Float'Value (Getenv (Name, Float'Image (Default)));
   end Getenv;

   function Getenv (Name : in String; Default : in Long_Integer) return Long_Integer is
   begin
      return Long_Integer'Value (Getenv (Name, Long_Integer'Image (Default)));
   end Getenv;

   function Getenv (Name : in String; Default : in Long_Float) return Long_Float is
   begin
      return Long_Float'Value (Getenv (Name, Long_Float'Image (Default)));
   end Getenv;

   function Getenv (Name : in String; Default : in Boolean) return Boolean is
   begin
      return Boolean'Value (Getenv (Name, Boolean'Image (Default)));
   end Getenv;

   procedure Setenv (Name : in String; Value : in String) is
   begin
      Ada.Environment_Variables.Set(Name, Value);
   end Setenv;

   procedure Clear (Name : in String) is
   begin
      Ada.Environment_Variables.Clear(Name);
   end Clear;

   -- String to list of integers
   function S_to_I (S : in string) return Ns is
     use Ada.Integer_Text_IO;
     List : Ns(1..100); -- magic number
     N, L : Integer := 0;
     Index : Integer := 1;
   begin
     loop
        get(S(L+1..S'Last), N, L);
        List(Index) := N;
        Index := Index+1;
     end loop;
   exception
     when others  =>
        return List(1..Index-1);
   end S_to_I;


   function S_to_LF (S : in string) return Fs is
     use Ada.Long_Float_Text_IO;
     List : Fs(1..100); -- magic number
     N : Long_Float;
     L : Integer := 0;
     Index : Integer := 1;
   begin
     loop
        get(S(L+1..S'Last), N, L);
        List(Index) := N;
        Index := Index+1;
     end loop;
   exception
     when others  =>
        return List(1..Index-1);
   end S_to_LF;


   function Command_Line_Option_Exists(Option : in String) return Boolean is
   begin
      for I in 1 .. Ada.Command_Line.Argument_Count loop
         if Ada.Command_Line.Argument(I) = Option then
            Text_IO.Put_Line("EXITING");
            return True;
         end if;
      end loop;
      return False;
   end Command_Line_Option_Exists;


   procedure Read_Response_File is
      FN : constant String := Ada.Directories.Simple_Name (Ada.Command_Line.Command_Name) & ".resp";
      use type Ada.Command_Line.Response_File.String_Access;
   begin
      Text_IO.Put_Line("RESPONSE FILE");
      declare
         L : Ada.Command_Line.Response_File.Argument_List :=
           Ada.Command_Line.Response_File.Arguments_From(FN);
      begin
         --Ada.Text_IO.Put_Line(L(2*I-1).all & " " & L(2*I).all);
         for I in L'First .. L'Last/2 loop
            for EV in Options loop
               if EV'Image = L(2*I-1).all then
                  OL(EV).Name :=  L(2*I-1);
                  OL(EV).Value :=  L(2*I);
               end if;
            end loop;
         end loop;
      end;
      for EV in Options loop
         if OL(EV).Name = null then
            Text_IO.Put_Line(EV'Image & " not set, use default");
         else
            Text_IO.Put_Line(OL(EV).Name.all & " " & OL(EV).Value.all);
         end if;
      end loop;

   exception
      when Ada.Command_Line.Response_File.File_Does_Not_Exist =>
         Text_IO.Put_Line(FN & " not found");
   end;

begin

   Read_Response_File;

end GEM;
