clc;
clear;

for lamb = [6000, 7000, 8000, 9000, 10000, 20000, 30000]
    genLambdaMUL(lamb);
end

function [] = genLambdaMUL(lamb)
    disp("MUL_lamb:" );
    disp(lamb)
    width = [8 8];
    g = width(1);
    h = width(2);
    l = 2;
    sign = 0;
    dN = 10;
    step = 100;
    
    if sign == 0
        mat_file = "unsigned_exchange_";
    else
        mat_file = "signed_exchange_";
    end
    mat_file = mat_file + string(g) + "x" + string(h) + "_l" +string(l);
    load(mat_file);
    
    objectFun = @(x) [1 x] * objectMat * [1 x]';
    
    % Z = ( 3*h - 2) * l;
    Z = ( 6*h - 8) * l;
    % upper bounds and lower bounds
    Lb = zeros(1, Z);
    Ub = ones(1, Z);
    % numbers listed in the IntCon takes integer values
    IntCon = 1:Z;
    
    fun = @(x) objectFun(x) + x*ones(Z, 1)*lamb;
    results = [];
    counter = 0;
    % options = optimoptions('ga','UseParallel', true, 'UseVectorized', false);
    % p = parpool(53);
    while counter <= 301
        disp(counter);
        % fval: the value of the fitness function at x.
        [x, fval] = ga(fun, Z, [], [], [], [], Lb, Ub, [], IntCon);
        results = [results; [sum(x) fval x]];
        results = unique(results, 'rows', 'stable');
        results = sortrows(results, 2);
        results = results(1:min(size(results,1), dN), :);
        if mod(counter, step) == 0
            for jdx = 1:min(size(results, 1), dN)
                genVerilog(g, h, l, lamb, jdx, results(jdx, :), dN, sign);
                genCpp(g, h, l, lamb, jdx, results(jdx, :), dN, sign);
            end
        end
        counter = counter + 1;
    end
end

function [] = genVerilog(g, h, l, lamb, file_number, result, dN, sign)
    verilog = "";
    verilog = verilog + "// terms: " + string(result(1)) + newline;
    verilog = verilog + "// fval:  " + string(result(2)) + newline + newline;

    if sign == 0
        mul_name = "unsigned_exchange_";
    else
        mul_name = "signed_exchange_";
    end

    mul_name = mul_name + string(g) + "x" + string(h) + "_l" + string(l) + "_lamb" + string(lamb) + "_" + num2str(file_number-1, '%0' + string( max(ceil(log10(dN)), 1) ) + "d");
    verilog = verilog + "module " + mul_name + " (" + newline;
    verilog = verilog + char(9) + "input [" + string( g-1 ) + ":0] x," + newline;
    verilog = verilog + char(9) + "input [" + string( h-1 ) + ":0] y," + newline;
    verilog = verilog + char(9) + "output [" + string( g+h-1 ) + ":0] z" + newline + ");" + newline + newline;

    for idx = 1 : g
        if sign == 0
            verilog = verilog + "wire [" + string( h-1 ) + ":0] part" + string(idx) + " =  y & {" + string(h) + "{x[" + string( idx-1 ) + "]}};" + newline;
        else % Baugh-Wooley structure
            verilog = verilog + "wire [" + string( h-1 ) + ":0] tmp_part" + string(idx) + " =  y & {" + string(h) + "{x[" + string( idx-1 ) + "]}};" + newline;
            verilog = verilog + "wire [" + string( h-1 ) + ":0] part" + string(idx) + ";" + newline;
            %%%%%%%%
            verilog = verilog + "assign part" + string(idx) + "[" + string( h-2 ) + ":0] = ";
            if idx == g
                verilog = verilog + "~";
            end
            verilog = verilog + "tmp_part" + string(idx) + "[" + string( h-2 ) + ":0];" + newline;
            %%%%%%%%
            verilog = verilog + "assign part" + string(idx) + "[" + string( h-1 ) + ":" + string( h-1 ) + "] = ";
            if idx ~= g
                verilog = verilog + "~";
            end
            verilog = verilog + "tmp_part" + string(idx) + "[" + string( h-1 ) + ":" + string( h-1 ) + "];";
            %%%%%%%%
            verilog = verilog + newline + newline;
        end
    end
    
    contents = get_contents(h, result);
    
    verilog = verilog + newline;
    [~, new_part_number] = mode(contents( : , 4));
    for idx = 1 : new_part_number
        new_part_wdith = max(contents( : , 4));
        verilog = verilog + "wire [" + string( new_part_wdith ) + ":0] new_part" + string(idx) + ";" + newline;
        % wire [12:0] new_part1;
        for jdx = 0 : new_part_wdith
            verilog = verilog + "assign new_part" + string(idx) + "[" + string(jdx) + "] = ";
            % Lia means true of false Loc represents the exact position
            [Lia, Loc] = ismember(jdx, contents( : , 4));
            if Lia == 0
                verilog = verilog + "0;";
            else
                row_number = mode(Loc, size(contents, 1));
                part_index = contents(row_number, 1);
                width_index = contents(row_number, 2);
                operation = contents(row_number, 3);
                % ddx and cdx!
                % verilog = verilog + "part" + string(part_index) + "[" + string(width_index) + "]";
                if operation == 0
                    verilog = verilog + "part" + string(part_index) + "[" + string(width_index) + "]";
                    verilog = verilog + ";";
                elseif operation == 1 % P1 AND P3
                    % verilog = verilog + " & part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "];";
                    verilog = verilog + "part" + string(part_index) + "[" + string(width_index) + "]";
                    verilog = verilog + " & part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "];";
                elseif operation == 2 % P1 OR P3
                    verilog = verilog + "part" + string(part_index) + "[" + string(width_index) + "]";
                    verilog = verilog + " | part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "];";
                elseif operation == 3 % P1 XOR P3
                    verilog = verilog + "part" + string(part_index) + "[" + string(width_index) + "]";
                    verilog = verilog + " ^ part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "];";

                elseif operation == 4 % P2 AND P4
                    % verilog = verilog + " & part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "];";
                    verilog = verilog + "part" + string(part_index) + "[" + string(width_index-1) + "]";
                    verilog = verilog + " & part" + string( part_index+1 ) + "[" + string( width_index-2 ) + "];";
                elseif operation == 5 % P2 OR P4
                    verilog = verilog + "part" + string(part_index) + "[" + string(width_index-1) + "]";
                    verilog = verilog + " | part" + string( part_index+1 ) + "[" + string( width_index-2 ) + "];";
                elseif operation == 6 % P2 XOR P4
                    verilog = verilog + "part" + string(part_index) + "[" + string(width_index-1) + "]";
                    verilog = verilog + " ^ part" + string( part_index+1 ) + "[" + string( width_index-2 ) + "];";

                elseif operation == 7 % P1 AND P4
                    % verilog = verilog + " & part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "];";
                    verilog = verilog + "part" + string(part_index) + "[" + string(width_index) + "]";
                    verilog = verilog + " & part" + string( part_index+1 ) + "[" + string( width_index-2 ) + "];";
                elseif operation == 8 % P1 OR P4
                    verilog = verilog + "part" + string(part_index) + "[" + string(width_index) + "]";
                    verilog = verilog + " | part" + string( part_index+1 ) + "[" + string( width_index-2 ) + "];";
                elseif operation == 9 % P1 XOR P4
                    verilog = verilog + "part" + string(part_index) + "[" + string(width_index) + "]";
                    verilog = verilog + " ^ part" + string( part_index+1 ) + "[" + string( width_index-2 ) + "];";

                elseif operation == 10 % P2 AND P3
                    % verilog = verilog + " & part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "];";
                    verilog = verilog + "part" + string(part_index) + "[" + string(width_index-1) + "]";
                    verilog = verilog + " & part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "];";
                elseif operation == 11 % P2 OR P3
                    verilog = verilog + "part" + string(part_index) + "[" + string(width_index-1) + "]";
                    verilog = verilog + " | part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "];";
                elseif operation == 12 % P2 XOR P3
                    verilog = verilog + "part" + string(part_index) + "[" + string(width_index-1) + "]";
                    verilog = verilog + " ^ part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "];";
                % else 
                %     disp("debugging: operation ==")
                %     disp(operation)
                end
                contents(row_number, : ) = [];
            end
            verilog = verilog + newline;
        end
        verilog = verilog + newline;
    end

    % assign z
    if sign == 0
        if l < g
            verilog = verilog + "wire [" + string( h+g-l-1 ) + ":0] tmp_z = y*x[" + string( g-1 ) + ":" + string(l) + "];" + newline + newline;
            verilog = verilog + "assign z = {tmp_z, " + string(l) + "'d 0} +";
        else
            verilog = verilog + "assign z =";
        end
    else
        verilog = verilog + "assign z =";
        for idx = l+1 : g
                verilog = verilog + " {part" + string(idx) + ", " + string(idx-1) + "'d 0} +";
        end
        verilog = verilog + " {1'd 1, " + string(g-2) + "'d 0, 1'd 1, " + string(h) + "'d 0} +"; % two constants '1' in Baugh-Wooley multiplier
    end


    for idx = 1 : new_part_number
        if idx == 1
            verilog = verilog + " new_part" + string(idx);
        else
            verilog = verilog + " + new_part" + string(idx);
        end
    end
    verilog = verilog + ";" + newline + newline + "endmodule";
    
    fid=fopen(mul_name + ".v", "wt");
    fprintf(fid, '%s\n', verilog);
    fclose(fid);

end

function [] = genCpp(g, h, l, lamb, file_number, result, dN, sign)
    cpp = "";
    cpp = cpp + "// terms: " + string(result(1)) + newline;
    cpp = cpp + "// fval:  " + string(result(2)) + newline + newline;
    cpp = cpp + "#include <bitset>" + newline;
    cpp = cpp + "#include <cmath>" + newline + newline;
    cpp = cpp + "using namespace std;" + newline + newline;

    if sign == 0
        type = "uint";
        mul_name = "unsigned_exchange_";
    else
        type = "int";
        mul_name = "signed_exchange_";
    end

    mul_name = mul_name + string(g) + "x" + string(h) + "_l" + string(l) + "_lamb" + string(lamb) + "_" + num2str(file_number-1, '%0' + string( max(ceil(log10(dN)), 1) ) + "d");

    x_width = 2^( ceil( log2( g ) ) );
    y_width = 2^( ceil( log2( h ) ) );
    total_width = 2^( ceil( log2(x_width+y_width) ) );

    cpp = cpp + type + string(total_width) + "_t " + mul_name + "(" + type + string(total_width) + "_t i, " + type + string(total_width) + "_t j) {" + newline + newline;
    cpp = cpp + "bitset<" + string( g ) + "> x = i;" + newline + newline;

    for idx = 1 : g
        cpp = cpp + "bitset<" + string( h ) + "> part" + string(idx) + " = j*x[" +  string(idx-1) + "];" + newline;
        if sign == 1
            cpp = cpp + "part" + string(idx) + "[" + string(h-1) + "].flip();" + newline + newline;
            if idx == g
                cpp = cpp + "part" + string(idx) + ".flip();" + newline;
            end
        end
    end

    cpp = cpp + newline;

    contents = get_contents(h, result);

    if sign == 0
        if l < g
            cpp = cpp + type + string(total_width) + "_t tmp_z = ";
            for idx = l : g-1
                if idx == l
                    cpp = cpp + " j*x[" + string(idx) + "]*pow(2, " + string(idx) + ")";
                else
                    cpp = cpp + " + j*x[" + string(idx) + "]*pow(2, " + string(idx) + ")";
                end
            end
            cpp = cpp + ";" + newline + newline;
            cpp = cpp + type + string(total_width) + "_t z = tmp_z +";
        else
            cpp = cpp + type + string(total_width) + "_t z =";
        end
    else
        cpp = cpp + type + string(total_width) + "_t z =";
        for idx = l+1 : g
            for jdx = 0 : h-1
                if ( idx == l+1 ) && ( jdx == 0 )
                    cpp = cpp + " part" + string(idx) + "[" + string(jdx) + "]*pow(2, " + string( idx-1 + jdx ) + ")";
                else cpp = cpp + " + part" + string(idx) + "[" + string(jdx) + "]*pow(2, " + string( idx-1 + jdx ) + ")";
                end
            end
        end
        cpp = cpp + " - pow(2, " + string( g+h-1 ) + ") + pow(2, " + string(h) + ") +"; % two constants '1' Baugh-Wooley multiplier
    end

    for row_number = 1 : size(contents, 1)
        part_index = contents(row_number, 1);
        width_index = contents(row_number, 2);
        operation = contents(row_number, 3);
        power = contents(row_number, 4);
        
        if row_number == 1
            cpp = cpp + " (part";
        else
            cpp = cpp + " + (part";
        end
        
        if operation == 0
            cpp = cpp + string(part_index) + "[" + string(width_index) + "]";
            cpp = cpp + "*pow(2, " + string(power) + "))";
        elseif operation == 1 % P1 AND P3
            cpp = cpp + string(part_index) + "[" + string(width_index) + "]";
            cpp = cpp + "&part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "])*pow(2, " + string(power) + ")";
        elseif operation == 2 % P1 OR P3
            cpp = cpp + string(part_index) + "[" + string(width_index) + "]";
            cpp = cpp + "|part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "])*pow(2, " + string(power) + ")";
        elseif operation == 3 % P1 XOR P3
            cpp = cpp + string(part_index) + "[" + string(width_index) + "]";
            cpp = cpp + "^part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "])*pow(2, " + string(power) + ")";

        elseif operation == 4 % P2 AND P4
            cpp = cpp + string(part_index) + "[" + string(width_index-1) + "]";
            cpp = cpp + "&part" + string( part_index+1 ) + "[" + string( width_index-2 ) + "])*pow(2, " + string(power) + ")";
        elseif operation == 5 % P2 OR P4
            cpp = cpp + string(part_index) + "[" + string(width_index-1) + "]";
            cpp = cpp + "|part" + string( part_index+1 ) + "[" + string( width_index-2 ) + "])*pow(2, " + string(power) + ")";
        elseif operation == 6 % P2 XOR P4
            cpp = cpp + string(part_index) + "[" + string(width_index-1) + "]";
            cpp = cpp + "^part" + string( part_index+1 ) + "[" + string( width_index-2 ) + "])*pow(2, " + string(power) + ")";

        elseif operation == 7 % P1 AND P4
            cpp = cpp + string(part_index) + "[" + string(width_index) + "]";
            cpp = cpp + "&part" + string( part_index+1 ) + "[" + string( width_index-2 ) + "])*pow(2, " + string(power) + ")";
        elseif operation == 8 % P1 OR P4
            cpp = cpp + string(part_index) + "[" + string(width_index) + "]";
            cpp = cpp + "|part" + string( part_index+1 ) + "[" + string( width_index-2 ) + "])*pow(2, " + string(power) + ")";
        elseif operation == 9 % P1 XOR P4
            cpp = cpp + string(part_index) + "[" + string(width_index) + "]";
            cpp = cpp + "^part" + string( part_index+1 ) + "[" + string( width_index-2 ) + "])*pow(2, " + string(power) + ")";

        elseif operation == 10 % P2 AND P3
            cpp = cpp + string(part_index) + "[" + string(width_index-1) + "]";
            cpp = cpp + "&part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "])*pow(2, " + string(power) + ")";
        elseif operation == 11 % P2 OR P3
            cpp = cpp + string(part_index) + "[" + string(width_index-1) + "]";
            cpp = cpp + "|part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "])*pow(2, " + string(power) + ")";
        elseif operation == 12 % P2 XOR P3
            cpp = cpp + string(part_index) + "[" + string(width_index-1) + "]";
            cpp = cpp + "^part" + string( part_index+1 ) + "[" + string( width_index-1 ) + "])*pow(2, " + string(power) + ")";
        end
    end

    cpp = cpp + ";";
    cpp = cpp + newline + newline + "return z;" + newline + newline;
    cpp = cpp + "}";
    
    fid = fopen(mul_name + ".h", "wt");
    fprintf(fid, '%s\n', cpp);
    fclose(fid);
end


function contents = get_contents(h, result)
    terms = 0;
    if mod(h,2) == 0
        terms = 12*h-16;
    else
        terms = 12*h-10;
    end
    % find the 2rd demensions (column) of the result
    x = result( 3 : size(result, 2) );
    contents = [];
    for jdx = 1 : size(x, 2)
        % if chosen
        if x(jdx) == 1
            % section A
            if mod(jdx-1, terms) == 0
                part_index = 2*( floor( (jdx-1) / (terms) ) ) + 1;
                width_index = 0;
                operation = 0; % no operation
                power = part_index - 1;
            % section B
            elseif mod(jdx-1, terms) == (terms) - 1
                part_index = 2*( floor( (jdx-1) / (terms) ) ) + 1 + 1;
                width_index = h - 1;
                operation = 0; % no operation
                power = part_index - 1 + h - 1;
            % section C
            elseif mod(h,2) == 0 && mod(jdx-1, terms) <= terms-2 && mod(jdx-1, terms) >= terms -2 - 5
                part_index = 2*( floor( (jdx-1) / (terms) ) ) + 1;
                width_index = h - 1;
                % pos' = mod(jdx-1, 12*h-16) - 1 -12*(h - 2)
                operation = floor( (mod(jdx-1, terms) - 1 -12*(h - 2)) / 2 ) + 1;
                % power = bp + mod(pos', 2)
                power = part_index - 1 + h - 1 + mod( (mod(jdx-1, terms) - 1 -12*(h - 2)), 2 );
            else
                % rowj = floor( (jdx-1) / (12*h-16) );
                % part_index = 2 * rowj +1;
                part_index = 2*( floor( (jdx-1) / (terms) ) ) + 1;
                % colj = mod(jdx-1, 12*h-16);
                % width_index = ddx - 1;
                width_index = 2 * floor( ( mod(jdx-1, terms) - 1 ) / 24 ) + 2;
                % pos = mod(colj - 1, 24) ;
                % 1: AND, 2: OR, 3: XOR
                operation = 3 * floor(mod( mod(jdx-1, terms) - 1, 24) / 6) + floor(mod(mod( mod(jdx-1, terms) - 1, 24),6) / 2) + 1;
                % bp(basic power) = width_index + (part_index - 1)
                % bias is for correction
                bias = [1, 0, 0, 1];
                power = width_index + (part_index - 1) -1 + mod(mod(mod(jdx-1, terms) - 1, 24), 2) + bias(floor(mod(mod(jdx-1, terms) - 1, 24) / 6) + 1);
            end
            % for debugging
            % if operation >= 13
            %     disp("debug: operation == ")
            %     disp(operation)
            %     disp("debug: jdx == ")
            %     disp(jdx)
            %     disp("debug: colj == ")
            %     disp( mod(jdx-1, 12*h-16))                
            %     disp("debug: pos' == ")
            %     disp(mod(jdx-1, 12*h-16) - 1 -12*(h - 2))
            %     disp("debug: operation == ")
            %     disp( floor( mod(jdx-1, 12*h-16) - 1 -12*(h - 2) / 2 ) + 1)
            %     disp("-------------------------")
            % end
            contents = [contents ; [part_index, width_index, operation, power]];
        end
    end
end
