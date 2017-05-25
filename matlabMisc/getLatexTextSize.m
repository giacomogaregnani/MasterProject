function sizeinpt = getLatexTextSize(sizename, varargin)

    % Command             10pt    11pt    12pt
    % \tiny               5       6       6
    % \scriptsize         7       8       8
    % \footnotesize       8       9       10
    % \small              9       10      10.95
    % \normalsize         10      10.95   12
    % \large              12      12      14.4
    % \Large              14.4    14.4    17.28
    % \LARGE              17.28   17.28   20.74
    % \huge               20.74   20.74   24.88
    % \Huge               24.88   24.88   24.88

    output = [  5       6       6        ;
                7       8       8        ;
                8       9       10       ;
                9       10      10.95    ;
                10      10.95   12       ;
                12      12      14.4     ;
                 14.4    14.4    17.28   ;
                17.28   17.28   20.74    ;
                20.74   20.74   24.88    ;
                24.88   24.88   24.88       ];


    p = inputParser;
    addParamValue(p,'docclass', '10pt');
    addParamValue(p,'enhanced', 1);
    parse(p,varargin{:});
    
    docclass = p.Results.docclass;
    enhanced = p.Results.enhanced;
    
    switch docclass
        case '10pt'
            ind = 1;
        case '11pt'
            ind = 2;
        case '12pt'
            ind = 3;
    end
    
    switch sizename
        case 'tiny';            sizeinpt = output(1,ind);
        case 'scriptsize';      sizeinpt = output(2,ind);
        case 'footnotesize';    sizeinpt = output(3,ind);
        case 'small';           sizeinpt = output(4,ind);
        case 'normalsize';      sizeinpt = output(5,ind);
        case 'large';           sizeinpt = output(6,ind);
        case 'Large';           sizeinpt = output(7,ind);
        case 'LARGE';           sizeinpt = output(8,ind);
        case 'huge';            sizeinpt = output(9,ind);
        case 'Huge';            sizeinpt = output(10,ind);
        otherwise;              error('invalid input');
    end
            
   sizeinpt = enhanced*sizeinpt;
    
end

