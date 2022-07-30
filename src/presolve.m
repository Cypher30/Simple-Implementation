function [A1,b1,Ix,valx,flag] = presolve(A,b)
[~,R,P] = qr(A',0);
i = length(P);
while abs(R(i,i)) < power(10,-10)
    i = i - 1;
end
A = A(P(1,1:i),:);
b = b(P(1,1:i),:);
[row,col] = find(A');
Ix = [];
Irow = [];
valx = [];
cflag = 0;
while 1
    flag_row = 1;
    key_x = length(Ix) + 1;
    key_row = 1;
    count = 0;
    countd = 1;
    for i = 1:length(col)
        if col(i) == flag_row
            if ~ismember(row(i),Ix)
                count = count + 1;
            else
                countd = countd + 1;
            end
        else
            if count == 1 && ~ismember(row(i-countd),Ix)
                Ix(key_x,1) = row(i - countd);
                valx(key_x,1) = b(col(i - countd),1) / A(col(i - countd),row(i - countd));
                if valx(key_x,1) < 0
                    fprintf("The given LP is infeasible!\n");
                    flag = 1;
                    A1 = 0;
                    b1 = 0;
                    Ix = 0;
                    valx = 0;
                    return;
                end
                if ismember(flag_row,Irow)
                    Irow(Irow == flag_row) = [];
                end
                b = b - valx(key_x,1) * A(:,row(i - countd));
                key_x = key_x + 1;
                count = 1;
                countd = 1;
                flag_row = flag_row + 1;
                cflag = 1;
                continue;
            else
                if ~ismember(flag_row,Irow) && count > 1
                    Irow(key_row,1) = flag_row;
                    key_row = key_row + 1;
                end
                if ~ismember(row(i),Ix)
                    count = 1;
                else
                    count = 0;
                end
                countd = 1;
                flag_row = flag_row + 1;
            end
        end
    end
    if count == 1 && ~ismember(row(i + 1 - countd),Ix)
        Ix(key_x,1) = row(i + 1 - countd);
        valx(key_x,1) = b(col(i + 1 - countd),1) / A(col(i + 1 - countd),row(i + 1 - countd));
        if valx(key_x,1) < 0
            fprintf("The given LP is infeasible!\n");
            flag = 1;
            A1 = 0;
            b1 = 0;
            Ix = 0;
            valx = 0;
            return;
        end
        b = b - valx(key_x,1) * A(:,row(i));
        if ismember(flag_row,Irow)
            Irow(Irow == flag_row) = [];
        end
        cflag = 1;
    else
        if ~ismember(flag_row,Irow) && count > 1
            Irow(key_row,1) = flag_row;
        end
    end
    if cflag == 0
        break;
    end
    cflag = 0;
end
if isempty(Ix)
    A1 = A;
    b1 = b;
    Ix = [];
    valx = [];
    flag = 0;
    return;
end
tempA = A(Irow,:);
b1 = b(Irow,:);
S = size(tempA);
flag1 = 1;
for i = 1:S(1,2)
    if ismember(i,Ix)
        continue;
    else
        A1(:,flag1) = tempA(:,i);
        flag1 = flag1 + 1;
    end
end
flag = 0;
end