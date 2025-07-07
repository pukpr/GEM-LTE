with GEM.LTE.Primitives;
function GEM.Mix_Regression (Target : in LTE.Primitives.Data_Pairs;
                             Start, Finish : in Long_Float;
                             First_Index, Last_Index : in Integer := 0;
                             Exclude : in Boolean := False;
                             Display_Results : in Boolean := False) return LTE.Primitives.Data_Pairs; -- Gem.LTE.Period_Set;

