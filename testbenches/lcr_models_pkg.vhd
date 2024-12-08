LIBRARY ieee  ; 
    USE ieee.NUMERIC_STD.all  ; 
    USE ieee.std_logic_1164.all  ; 
    use ieee.math_real.all;

package lcr_models_pkg is

------------------------------------------------------
    type lcr_model_3ph_record is record
        states : real_vector(0 to 5);
        un : real;
    end record;

    function init_lcr_model return lcr_model_3ph_record ;

    function get_capacitor_voltage(lcr : lcr_model_3ph_record) 
        return real_vector;
    function get_inductor_current(lcr : lcr_model_3ph_record) 
        return real_vector;
------------------------------------------------------
    function deriv_lcr (
        states   : real_vector
        ; i_load : real_vector
        ; uin    : real_vector
        ; l      : real_vector
        ; c      : real_vector
        ; r      : real_vector)
        return lcr_model_3ph_record;

    function deriv_lcr (
        states   : real_vector
        ; i_load : real_vector
        ; uin    : real_vector
        ; l      : real_vector
        ; c      : real_vector
        ; r      : real_vector) 
        return real_vector ;

------------------------------------------------------
    function get_neutral_voltage(ul : real_vector; l : real_vector) return real;
------------------------------------------------------

end package;

package body lcr_models_pkg is

    use work.real_vector_pkg.all;

------------------------------------------------------
    function init_lcr_model return lcr_model_3ph_record 
    is
    begin
        return ((others => 0.0), 0.0);
    end init_lcr_model;
------------------------------------------------------
    function get_capacitor_voltage(lcr : lcr_model_3ph_record) 
        return real_vector
    is
        variable retval : real_vector(0 to 2);
    begin
        retval := lcr.states(3 to 5) - lcr.un;
        return retval;
    end get_capacitor_voltage;
------------------------------------------------------
    function get_inductor_current(lcr : lcr_model_3ph_record) 
        return real_vector
    is
        variable retval : real_vector(0 to 2);
    begin
        retval := lcr.states(0 to 2);
        return retval;
    end get_inductor_current;
------------------------------------------------------
    function get_neutral_voltage(ul : real_vector; l : real_vector) return real
    is
        constant div : real                := 1.0/(l(1)*l(2) + l(1)*l(3) + l(2)*l(3));
        constant a   : real_vector(1 to 3) := (l(2)*l(3), l(1)*l(3), l(1)*l(2));
    begin
        return (a(1)*ul(1) + a(2)*ul(2) + a(3)*ul(3))*div;
    end function;

------------------------------------------------------
    function calculate_lcr (
        states   : real_vector
        ; i_load : real
        ; uin    : real
        ; l      : real
        ; c      : real
        ; r      : real) 
        return real_vector is

        variable retval : real_vector(0 to 1);

        variable ul : real;
        alias il is states(0);
        alias uc is states(1);

        variable dil : real;
        variable duc : real;

    begin
        ul := uin - uc - il * r;
        dil := ul/l;
        duc := (il - i_load) / c;

        retval := (dil, duc);

        return retval;
    end calculate_lcr;
------------------------------------------------------
    function deriv_lcr (
        states   : real_vector
        ; i_load : real_vector
        ; uin    : real_vector
        ; l      : real_vector
        ; c      : real_vector
        ; r      : real_vector) 
        return lcr_model_3ph_record 
    is

        variable retval : lcr_model_3ph_record;

        variable ul : real_vector(1 to 3) := (0.0 , 0.0 , 0.0);
        alias il    : real_vector(1 to 3) is states(0 to 2);
        alias uc    : real_vector(1 to 3) is states(3 to 5);

        variable un  : real := 0.0;

        variable dil : real_vector(1 to 3);
        variable duc : real_vector(1 to 3);

    begin
        ul(1) := uin(1) - uc(1) - il(1) * r(1);
        ul(2) := uin(2) - uc(2) - il(2) * r(2);
        ul(3) := uin(3) - uc(3) - il(3) * r(3);
        un := get_neutral_voltage(ul, l);

        dil(1) := (ul(1)-un)/l(1);
        dil(2) := (ul(2)-un)/l(2);
        dil(3) := (ul(3)-un)/l(3);

        duc(1) := (il(1) - i_load(0))               / c(1);
        duc(2) := (il(2) - i_load(1))               / c(2);
        duc(3) := (il(3) - (i_load(1) - i_load(0))) / c(3);

        retval := ( ( dil(1), dil(2), dil(3), duc(1), duc(2), duc(3) )
                    , un);

        return retval;

    end deriv_lcr;

    function deriv_lcr (
        states   : real_vector
        ; i_load : real_vector
        ; uin    : real_vector
        ; l      : real_vector
        ; c      : real_vector
        ; r      : real_vector) 
        return real_vector 
    is
        variable retval : lcr_model_3ph_record;
    begin
        retval := deriv_lcr(states, i_load, uin, l, c, r);
        return retval.states;

    end deriv_lcr;






end package body;
-----------------
