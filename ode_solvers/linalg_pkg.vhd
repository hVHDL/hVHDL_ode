LIBRARY ieee  ;
    USE ieee.NUMERIC_STD.all  ;
    USE ieee.std_logic_1164.all  ;
    use ieee.math_real.all;

-- Minimal dense linear algebra used by the implicit ODE solvers.
package linalg_pkg is

    type real_matrix is array(natural range <>, natural range <>) of real;

    -- Solves A*x = b by Gaussian elimination with partial pivoting.
    -- A must be square and indexed (b'range, b'range).
    function solve(A : real_matrix; b : real_vector) return real_vector;

end package linalg_pkg;

package body linalg_pkg is

    function solve(A : real_matrix; b : real_vector) return real_vector is
        constant lo : integer := b'low;
        constant hi : integer := b'high;
        variable m  : real_matrix(lo to hi, lo to hi);
        variable r  : real_vector(lo to hi);
        variable x  : real_vector(lo to hi) := (others => 0.0);
        variable pivot_row : integer;
        variable maxabs    : real;
        variable factor    : real;
        variable tmp       : real;
    begin
        for i in lo to hi loop
            r(i) := b(i);
            for j in lo to hi loop
                m(i, j) := A(i, j);
            end loop;
        end loop;

        -- forward elimination
        for k in lo to hi loop
            -- partial pivot
            pivot_row := k;
            maxabs    := abs(m(k, k));
            for row in k + 1 to hi loop
                if abs(m(row, k)) > maxabs then
                    maxabs    := abs(m(row, k));
                    pivot_row := row;
                end if;
            end loop;
            if pivot_row /= k then
                for col in lo to hi loop
                    tmp                := m(k, col);
                    m(k, col)          := m(pivot_row, col);
                    m(pivot_row, col)  := tmp;
                end loop;
                tmp          := r(k);
                r(k)         := r(pivot_row);
                r(pivot_row) := tmp;
            end if;

            for row in k + 1 to hi loop
                factor := m(row, k) / m(k, k);
                for col in k to hi loop
                    m(row, col) := m(row, col) - factor * m(k, col);
                end loop;
                r(row) := r(row) - factor * r(k);
            end loop;
        end loop;

        -- back substitution
        for k in hi downto lo loop
            tmp := r(k);
            for col in k + 1 to hi loop
                tmp := tmp - m(k, col) * x(col);
            end loop;
            x(k) := tmp / m(k, k);
        end loop;

        return x;
    end function solve;

end package body linalg_pkg;
