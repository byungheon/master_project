function joint=stateTojoint(state)
    if(size(state) == [3 1])
        joint       = state;
        joint(3)    = joint(3) - joint(1) + joint(2);
    else
        error('state converting error');
    end
end