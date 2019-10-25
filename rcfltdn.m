function new_signal = rcfltdn(raw_signal)
%RCFLTDN 此处显示有关此函数的摘要
%   此处显示详细说明
%     raw_signal = 1:20000;
%     b = rcosdesign(0.5, 4, 10);
%     new_signal = upfirdn(raw_signal, b, 1, 10);
    beta = 10;
    Nsym = 4;
    sampsPerSym = 10;
    H = [0.013427674865534088, 0.012735903231668802,...
        0.008453156912747329, 0.00043693057910991815, -0.010739957583904041,...
        -0.023736999882468295, -0.03651359231196992, -0.046520709600807369,...
        -0.051000786326600635, -0.047360864551823735, -0.03356918716383521,...
        -0.00851745396358352, 0.027708352033020543, 0.073695928108432215,...
        0.12674306972477839, 0.18306897632094107, 0.23818153519927915,...
        0.28735810559134217, 0.32617926034852995, 0.35104533721079589,...
        0.35960619065433957, 0.35104533721079589, 0.32617926034852995,...
        0.28735810559134217, 0.23818153519927915, 0.18306897632094107,...
        0.12674306972477839, 0.073695928108432215, 0.027708352033020543,...
        -0.00851745396358352, -0.03356918716383521, -0.047360864551823735,...
        -0.051000786326600635, -0.046520709600807369, -0.03651359231196992,...
        -0.023736999882468295, -0.010739957583904041, 0.00043693057910991815,...
        0.008453156912747329, 0.012735903231668802, 0.013427674865534088];
    StartIdx = [0, 1, 1, 1, 1, 1, 1, 1, 1, 1];
    
    N = length(raw_signal);
    new_signal = zeros(fix(N/sampsPerSym)+5, 1);
    InBuf = zeros(41, 1);
    
    outIdx = 1;
    inIdx = 1;
    inBufIdx = 1;
    for n = 1:N
        nModDFactor = mod(n-1, 10);
        k = StartIdx(nModDFactor+1);
        nModDFactor = 1;
        
        InBuf(inBufIdx) = raw_signal(inIdx);
        inIdx = inIdx + 1;
        
        while k < nModDFactor
            acc = 0.0;
%             coefPolyphaseOffset = b_obj->P1_PolyphaseSelector * 41;
            for i = inBufIdx:41
                acc = acc + H(i - inBufIdx + 1) * InBuf(i);
            end
            for i = 1:inBufIdx-1
                acc = acc + H(i - inBufIdx + 41 + 1) * InBuf(i);
            end

            new_signal(outIdx) = acc;
            outIdx = outIdx + 1;
            k = k + 1;
        end
        
        if inBufIdx == 1
            inBufIdx = 41;
        else
            inBufIdx = inBufIdx - 1;
        end
    end
end

