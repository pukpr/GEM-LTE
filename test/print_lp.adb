--  print_lp.adb - Print LP periods from GEM.LTE
--
--  Outputs the predefined LP (Long Period) tidal constituent periods
--  computed from Doodson arguments. These are the canonical periods
--  used by the JSON reader's period-matching logic.

with Ada.Text_IO;            use Ada.Text_IO;
with Ada.Long_Float_Text_IO; use Ada.Long_Float_Text_IO;
with GEM.LTE;

procedure Print_LP is
begin
   Put_Line ("LP_PERIODS = [");
   for I in GEM.LTE.LP'Range loop
      Put ("    ");
      Put (GEM.LTE.LP (I), Exp => 0, Aft => 11);
      Put_Line (",");
   end loop;
   Put_Line ("]");
end Print_LP;
