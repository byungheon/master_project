function state=jointTostate(joint)
    if(size(joint) == [3 1])
        state       = joint;
        state(3)    = state(3) + state(1) - state(2);
    else
        error('state converting error');
    end
end