function printQuadTriangleCase(q, caseNum)

fprintf('    case %d\n', caseNum);

fprintf('        q.Points =     [');
P = q.Points;
for i = 1:size(P,1)
    fprintf('%20.15f   %20.15f', P(i,1), P(i,2));
    if i < size(P,1)
        fprintf('\n                           ');
    end
end
fprintf('];\n');

fprintf('        q.Weights =    [');
W = q.Weights;
for i = 1:length(W)
    fprintf('%20.15f', W(i));
    if i < length(W)
        fprintf('\n                           ');
    end
end
fprintf('];\n');

end