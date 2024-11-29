LIBRARY ieee  ; 
    USE ieee.NUMERIC_STD.all  ; 
    USE ieee.std_logic_1164.all  ; 
    use ieee.math_real.all;

package lcr_models_pkg is

    impure function deriv_lcr (
        states : real_vector
        ; i_load : real_vector
        ; uin : real_vector
        ; c : real_vector
        ; l : real_vector 
        ; r : real_vector) 
        return real_vector;

end package;

package body lcr_models_pkg is

    impure function deriv_lcr (
        states : real_vector
        ; i_load : real_vector
        ; uin : real_vector
        ; c : real_vector
        ; l : real_vector
        ; r : real_vector) 
        return real_vector is

        variable retval : real_vector(0 to 5);

        variable ul : real_vector(1 to 3) := (0.0 , 0.0 , 0.0);
        alias il    : real_vector(1 to 3) is states(0 to 2);
        alias uc    : real_vector(1 to 3) is states(3 to 5);

        variable un  : real := 0.0;

        constant div : real                := 1.0/(l(1)*l(2) + l(1)*l(3) + l(2)*l(3));
        constant a   : real_vector(1 to 3) := (l(2)*l(3)/div, l(1)*l(3)/div, l(1)*l(2)/div);


        variable dil : real_vector(1 to 3);
        variable duc : real_vector(1 to 3);

    begin
        ul(1) := uin(1) - uc(1) - il(1) * r(1);
        ul(2) := uin(2) - uc(2) - il(2) * r(2);
        ul(3) := uin(3) - uc(3) - il(3) * r(3);
        un := a(1)*ul(1) + a(2)*ul(2) + a(3)*ul(3);

        dil(1) := (ul(1)-un)/l(1);
        dil(2) := (ul(2)-un)/l(2);
        dil(3) := (ul(3)-un)/l(3);

        duc(1) := (il(1) - i_load(0))/c(1);
        duc(2) := (il(2) - i_load(1))/c(2);
        duc(3) := (il(3) - (i_load(1) - i_load(0)))/c(3);

        retval := (dil(1), dil(2), dil(3), duc(1), duc(2), duc(3));

        return retval;

    end deriv_lcr;

end package body;
-----------------
