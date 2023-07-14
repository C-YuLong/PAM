

clear;
clc;

width = [8 8]; % width of a multiplier
g = width(1); % a multiplier g rows of partial products
h = width(2); % each row contains h partial products
l = 6; % select l rows to be compressed
sign = 0; % unsign or sign, 0: unsign, 1: sign

%%
 % Z compression term HERE! we assume that h is even!
if mod(h,2) == 0
    Z = l*(6*h-8);
else
    Z = l*(6*h-5);
end
tmpMat = zeros(Z+1, Z+1);
objectMat = zeros(Z+1, Z+1);


if sign == 0
    x_start = 0;
    x_end = int64(2^g-1);
    y_start = 0;
    y_end = int64(2^h-1);
else
    x_start = int64(-2^(g-1));
    x_end = int64(2^(g-1)-1);
    y_start = int64(-2^(h-1));
    y_end = int64(2^(h-1)-1);
end


% p = parpool(128);
tFor = tic; %start the timer
% calculating the partial product
for idx = x_start : x_end
%     tic;
    % display(idx);
    n = zeros(Z+1, 1);
    part_pro =  false(g ,h) ;
    for jdx = y_start : y_end
        % partial products
        % calculating the partial products
        for cdx = 1:g
            for ddx = 1:h
                part_pro(cdx, ddx) = bitget(jdx, ddx) * bitget(idx, cdx);
                if sign == 1
                    if ( cdx ~= g &&  ddx == h ) || ( cdx == g &&  ddx ~= h )
                        part_pro(cdx, ddx) = ~part_pro(cdx, ddx);
                    end
                end
            end
        end
% uncomressed partial products
n(1) = idx * jdx;
for cdx = l+1 : g
    for ddx = 1 : h
        % n(1) is the exact product of the idx and jdx
        % so we have to substract the uncompressed partial products
        % from the exact result
        n(1) = n(1) - part_pro(cdx, ddx) * 2^( (cdx-1) + (ddx-1) );
    end
end
if sign == 1
    n(1) = n(1) + 2^( g+h-1 ) - 2^h; % two constant '1' in Baugh-Wooley multiplier
end
% compression term candidates
counter = 2;
for cdx = 1 : 2 : l
    for ddx = 1 : h
    % to all the conditions, wheather odd or even
        if ddx == 1
            % the redundant term 1
            n(counter) = (-1) * part_pro(cdx, ddx) * 2^( (cdx-1) + (ddx-1) );
            counter = counter + 1
            ddx = ddx + 2
        elseif (ddx == h)
            % the redundant term 2
            n(counter) = (-1) * part_pro(cdx+1, ddx) * 2^( (cdx) + (ddx-1) );
            counter = counter + 1
        % if h is even, the last column should be considered separately
        elseif (mod(h,2) == 0) && (ddx == h)
            % AND
            n(counter) = (-1) * ( part_pro(cdx, h) & part_pro(cdx+1, h-1) ) * 2^( (cdx-1) + (h-1) );
            counter = counter + 1;
            n(counter) = (-1) * ( part_pro(cdx, h) & part_pro(cdx+1, h-1) ) * 2^( (cdx-1) + (h-1) + 1 );
            counter = counter + 1;
            % OR
            n(counter) = (-1) * ( part_pro(cdx, h) | part_pro(cdx+1, h-1) ) * 2^( (cdx-1) + (h-1) );
            counter = counter + 1;
            n(counter) = (-1) * ( part_pro(cdx, h) | part_pro(cdx+1, h-1) ) * 2^( (cdx-1) + (h-1) + 1 );
            counter = counter + 1;
            % XOR
            n(counter) = (-1) * xor( part_pro(cdx, h), part_pro(cdx+1, h-1) ) * 2^( (cdx-1) + (h-1) );
            counter = counter + 1;
            n(counter) = (-1) * xor( part_pro(cdx, h), part_pro(cdx+1, h-1) ) * 2^( (cdx-1) + (h-1) + 1 );
            counter = counter + 1;
        else
            % P1 AND P3
            n(counter) = (-1) * ( part_pro(cdx, ddx) & part_pro(cdx+1, ddx-1) ) * 2^( (cdx-1) + (ddx-1) );
            counter = counter + 1;
            n(counter) = (-1) * ( part_pro(cdx, ddx) & part_pro(cdx+1, ddx-1) ) * 2^( (cdx-1) + (ddx-1) + 1 );
            counter = counter + 1;
            % P1 OR P3
            n(counter) = (-1) * ( part_pro(cdx, ddx) | part_pro(cdx+1, ddx-1) ) * 2^( (cdx-1) + (ddx-1) );
            counter = counter + 1;
            n(counter) = (-1) * ( part_pro(cdx, ddx) | part_pro(cdx+1, ddx-1) ) * 2^( (cdx-1) + (ddx-1) + 1 );
            counter = counter + 1;
            % P1 XOR P3
            n(counter) = (-1) * xor( part_pro(cdx, ddx), part_pro(cdx+1, ddx-1) ) * 2^( (cdx-1) + (ddx-1) );
            counter = counter + 1;
            n(counter) = (-1) * xor( part_pro(cdx, ddx), part_pro(cdx+1, ddx-1) ) * 2^( (cdx-1) + (ddx-1) + 1 );
            counter = counter + 1;

            % P2 AND P4
            n(counter) = (-1) * ( part_pro(cdx, ddx-1) & part_pro(cdx+1, ddx-2) ) * 2^( (cdx-1) + (ddx-2) );
            counter = counter + 1;
            n(counter) = (-1) * ( part_pro(cdx, ddx-1) & part_pro(cdx+1, ddx-2) ) * 2^( (cdx-1) + (ddx-2) + 1 );
            counter = counter + 1;
            % P2 OR P4
            n(counter) = (-1) * ( part_pro(cdx, ddx-1) | part_pro(cdx+1, ddx-2) ) * 2^( (cdx-1) + (ddx-2) );
            counter = counter + 1;
            n(counter) = (-1) * ( part_pro(cdx, ddx-1) | part_pro(cdx+1, ddx-2) ) * 2^( (cdx-1) + (ddx-2) + 1 );
            counter = counter + 1;
            % P2 XOR P4
            n(counter) = (-1) * xor( part_pro(cdx, ddx-1), part_pro(cdx+1, ddx-2) ) * 2^( (cdx-1) + (ddx-2) );
            counter = counter + 1;
            n(counter) = (-1) * xor( part_pro(cdx, ddx-1), part_pro(cdx+1, ddx-2) ) * 2^( (cdx-1) + (ddx-2) + 1 );
            counter = counter + 1;

            % P1 AND P4
            n(counter) = (-1) * ( part_pro(cdx, ddx) & part_pro(cdx+1, ddx-2) ) * 2^( (cdx-1) + (ddx-2) );
            counter = counter + 1;
            n(counter) = (-1) * ( part_pro(cdx, ddx) & part_pro(cdx+1, ddx-2) ) * 2^( (cdx-1) + (ddx-2) + 1 );
            counter = counter + 1;
            % P1 OR P4
            n(counter) = (-1) * ( part_pro(cdx, ddx) | part_pro(cdx+1, ddx-2) ) * 2^( (cdx-1) + (ddx-2) );
            counter = counter + 1;
            n(counter) = (-1) * ( part_pro(cdx, ddx) | part_pro(cdx+1, ddx-2) ) * 2^( (cdx-1) + (ddx-2) + 1 );
            counter = counter + 1;
            % P1 XOR P4
            n(counter) = (-1) * xor( part_pro(cdx, ddx), part_pro(cdx+1, ddx-2) ) * 2^( (cdx-1) + (ddx-2) );
            counter = counter + 1;
            n(counter) = (-1) * xor( part_pro(cdx, ddx), part_pro(cdx+1, ddx-2) ) * 2^( (cdx-1) + (ddx-2) + 1 );
            counter = counter + 1;
            
            % P2 AND P3
            n(counter) = (-1) * ( part_pro(cdx, ddx-1) & part_pro(cdx+1, ddx-1) ) * 2^( (cdx-1) + (ddx-2) );
            counter = counter + 1;
            n(counter) = (-1) * ( part_pro(cdx, ddx-1) & part_pro(cdx+1, ddx-1) ) * 2^( (cdx-1) + (ddx-2) + 1 );
            counter = counter + 1;
            % P2 OR P3
            n(counter) = (-1) * ( part_pro(cdx, ddx-1) | part_pro(cdx+1, ddx-1) ) * 2^( (cdx-1) + (ddx-2) );
            counter = counter + 1;
            n(counter) = (-1) * ( part_pro(cdx, ddx-1) | part_pro(cdx+1, ddx-1) ) * 2^( (cdx-1) + (ddx-2) + 1 );
            counter = counter + 1;
            % P2 XOR P3
            n(counter) = (-1) * xor( part_pro(cdx, ddx-1), part_pro(cdx+1, ddx-1) ) * 2^( (cdx-1) + (ddx-2) );
            counter = counter + 1;
            n(counter) = (-1) * xor( part_pro(cdx, ddx-1), part_pro(cdx+1, ddx-1) ) * 2^( (cdx-1) + (ddx-2) + 1 );
            counter = counter + 1;
        end
    end
end

% %         n = n * ((gw(1, idx+1)/prow)^(0.5)) * ((gl(1, jdx+1)/prol)^(0.5)); % gw gl
% %         n = n * ((gw(1, idx + 2^(bit-1) + 1) / prow)^(0.5)); % FIR weight signal
        n = n / ( (2^g)^0.5 ) / ( (2^h)^0.5 ); % uniform distribution
        tmpMat = n * n.';
        objectMat = objectMat + tmpMat;
    end
    % toc;
end
display(toc(tFor));
% delete(p);

if sign == 0
    file_name = "unsigned_";
else
    file_name = "signed_";
end

file_name = file_name + string(g) + "x" + string(h) + "_l" +string(l);
save(file_name, 'objectMat');

% quit;
        