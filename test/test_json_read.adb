--  test_json_read.adb - Test JSON LPAP reading
--
--  Simple test to verify LPAP data is being read correctly from JSON

with Ada.Text_IO;            use Ada.Text_IO;
with Ada.Integer_Text_IO;    use Ada.Integer_Text_IO;
with Ada.Long_Float_Text_IO; use Ada.Long_Float_Text_IO;
with GEM.LTE.Primitives.Shared;
with GEM.LTE;

procedure Test_JSON_Read is
   Params : GEM.LTE.Primitives.Shared.Param_S (29, 11);
   Success : Boolean;
begin
   Success := GEM.LTE.Primitives.Shared.Read_JSON ("lt.exe.json", Params, True);
   
   if Success then
      Put_Line ("JSON read successful!");
      Put_Line ("LPAP data (first 5):");
      for I in 1 .. 5 loop
         Put ("  ");
         Put (I, Width => 2);
         Put (". Period=");
         Put (Params.A.LP (I), Exp => 0, Aft => 11);
         Put (" Amp=");
         Put (Params.B.LPAP (I).Amplitude, Exp => 0, Aft => 11);
         Put (" Phase=");
         Put (Params.B.LPAP (I).Phase, Exp => 0, Aft => 11);
         New_Line;
      end loop;
      Put_Line ("harm data (first 10):");
      for I in 1 .. 10 loop
         Put ("  ");
         Put (I, Width => 2);
         Put (". ");
         Put (Params.C (I), Width => 5);
         New_Line;
      end loop;
   else
      Put_Line ("JSON read failed!");
   end if;
end Test_JSON_Read;
