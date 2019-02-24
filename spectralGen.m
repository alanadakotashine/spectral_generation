function spectralGen(graphFile,outputFile,numOptSteps,numGraphs,optFlag)
%-------------------------------------------------------------
%Given an input graph,  
%generates graph with symmetric normalized laplacian
%whose eigenvalues are
%approximate to the target.  

%graphFile is the path to a gml file
%outputFile is the path to directory in which to place generated graphs
%numOptSteps is how many maximum steps to take in stiefel manifold
%optimization
%numRounds is how many graphs to generate
%resulting graphs will be output to outputFile in a csv accomanied by the
%spectrum of their symmetric normalized laplacian in a txt file
%-------------------------------------------------------------

%read in the graph
graph = read_gml(graphFile);
A = get_adjacency_matrix(graph);

A = A + A';
%reomve self loops
A = A - diag(diag(A));
n = size(A,1);
degreeSequence = sum(A);
targetNumEdges = sum(degreeSequence);
deg = diag(degreeSequence);
degSqrt = sqrt(deg);
invdegSqrt = diag(sqrt(diag(deg)).^-1);
L = deg - A;
Lap = invdegSqrt*L*invdegSqrt;
[trueX,Lambda] = eig(Lap);
targetSpec = diag(Lambda);
[targetSpec,order] = sort(targetSpec);

[fiedlerCut_truth,fiedlerCutComp_truth] = genFiedCut(trueX,order);
trueCutSize = min(size(fiedlerCut_truth,2),size(fiedlerCutComp_truth,2));


for numSucc = 1:numGraphs
    %file naming
    loc = strcat(outputFile,'/');
    %generate starting basis using the configuration model
    [startX,spectrumConfigModel,configGraph,raw] = configurationModelBenchmark(degreeSequence);
    %template matrix
    G = genTemplate(startX); 
    %iterate on template matrix if optimization flag is on by optimizing
    %over stiefel manifold to find better basis
    if(optFlag == 1)
        [X,G] = lineSearch(startX,@binary, @polarObj, @polarStep, @polarDer, @steepestDescentPolar, 100, .01, .0001,.01,numOptSteps);
    end
    %find fractional matrix while presrving degrees 
    [G] = genFractional(G);
    %drive down spectral gap by coordinating edge additions and deletions
    %according to fiedler cut
    [G] = driveDownSpecGap(G,targetSpec);
    %indpendently flip all remaining edges
    G = singleFlipAll(G);
    %get rid of all floats
    G(G<.01)=0;
    G(G>.9)=1;
    %record the result
    recordSpec(G,'specGen','specGenSpectrum',loc,numSucc);
    clear G
end

%%%%%%%%%%%Linear programs

    function [Y] = genTemplate(Q)
        LQ = Q*Lambda*Q';
        %n variables for square root degrees-x
        %n variables for result of LQx
        %b = LQx
        %minimize b st x > epsilon
        epsilon = .0001;
        AIn = [zeros(n),eye(n),-1*ones(n,1);zeros(n),-1*eye(n),-1*ones(n,1)];
        Aeq = [LQ,-1*eye(n),zeros(n,1)];
        lb = [ones(1,n)*epsilon,ones(1,n+1)*-1*Inf];
        ub = ones(1,2*n+1)*Inf;
        [x,fval,exitflag,output,lambda] = cplexlp([zeros(2*n,1);1],AIn,zeros(2*n,1),Aeq,zeros(n,1),lb,ub);
        ds_v = x(1:n);
        t=ds_v.*ds_v;
        ds=diag(ds_v);
        degCur = diag(t);
        tmp = ds*LQ*ds;
        Y = degCur - tmp;
        Y = real(Y);
        %rescale Y so that it has the correct number of edges
        scale = targetNumEdges/sum(sum(Y));
        Y=Y.*scale;
    end

    function[A] = genFractionalLp(Y,deq)
        %2q constraints to keep A_ij+E_ij between 0 and 1
        %q constraints to keep E_ij < G_ij
        %q constraints to keep -E_ij < G_ij
        %n equality constraints for degree or 2n inequality constraints for
        %degree
        %n degree constraints
        q = (n^2-n)/2;
        numVars=2*q;
        f = zeros(1,numVars);
        k=1;
        %penalties proportional to degrees
        for i = 1:n
            degi = deq(i);
            for j = i+1:n
                f(q+k) = 1/sqrt(degi*deq(j));
                k=k+1;
            end
        end
        Y=real(Y);
        Y = Y-diag(diag(Y));
        [AIn,b,Aeq,beq] = roundLinProgPrev(Y,deq);
        lb = [ones(1,q)*-Inf,zeros(1,q)];
        ub = [ones(1,numVars)*Inf];
        [x,fval,exitflag] = cplexlp(f',AIn,b',Aeq,beq',lb,ub);
        tri = triu(ones(n,n),1)';
        E = tri;
        E(tri==1)=x(1:q);
        E = E+E';
        A=Y+E;
    end

    function[AIn,b,Aeq,beq] = roundLinProgPrev(Y,degrees)
        %graph Y passed in has a zero diagonal
        %degrees past in is target degree sequence
        %one variable for each adjacency matrix entry corresponding to
        %a proper edge
        numEdges = (n^2-n)/2;
        %q variables for epsilon, q variables for Q, n variables for self-loops,
        %n variables degree change
        %2q constraints to keep E_ij < Q_ij and -E_ij < Q_ij, computes aboslute
        %difference Q_ij
        %2q constraints to keep Y_ij + E_ij between 0 and 1
        %n constraints equality for the degree
        sparseI = speye(numEdges);
        sparseNegI = -1*sparseI;
        sparseZero = sparse(numEdges,numEdges);
        AIn=[sparseI,sparseNegI;sparseNegI,sparseNegI;sparseI,sparseZero;sparseNegI,sparseZero];
        emptyTri = (triu(ones(n,n),1))';
        yFlat = Y(emptyTri==1);
        b = sparse([zeros(1,2*numEdges),1-yFlat',yFlat']);
        %computes difference between target and our current degree sequence
        beq=degrees-sum(Y);
        hack = [n];
        %there are n constraints, and 2*numEdges variables. The absolute
        %value difference variables corresponding to each edges are not included in
        %the equality constraints. C code generates equality constraints on
        %the difference variables and returns matrix.
        %this is appended to a zero matrix to form Aeq
        Aeq = [genEqualityConstraints(hack),sparse(n,numEdges)];
    end







   






    function[Y] = singleFlipAll(Y,nonBinPairs)
        rng('shuffle')
        ut = triu(Y,0);
        ut(ut<0)=0.0;
        ut(ut>1)=1.0;
        ut = binornd(1,ut);
        Y = ut'+ut;
    end




    



    function[vec,sizeCut,k] = findValidCut(minSize,normLap)
        [Xy,curSpec] = eig(normLap);
        curSpec = diag(curSpec);
        [useless,perm_1]=sort(curSpec);
        Xy=Xy(:,perm_1);
        k = 2;
        while true
            vec = Xy(:,k);
            [posC,negC] =genCut(vec);
            sizeCut = min(size(posC,2),size(negC,2));
            if sizeCut >= minSize
                break
            end
            k=k+1;
        end
        clear Xy
    end

    function[rG] = criticalCutRounding(rG,targetSpec)
        normLapY = computeNormalizedLaplacian(rG);
        [e2,lambda2] = genSecondCutAndEig(normLapY);
        count = 0;
        amountRemovedTotal = 0;
        [v,sizeCut,k] = findValidCut(.5*trueCutSize,normLapY);
        %disp('alter spec gap, should be negative otherwise it wont work');
        delta = targetSpec(k)-kthEig(normLapY,k);
        threshold = .5*abs(delta);
        numEdgesCrossing = compCrossEdges(v,rG);
        budgetParam = .25;
        if delta < 0
            budget = numEdgesCrossing*budgetParam;
        end
        numEdgesCrossing_fiedler = compCrossEdges(e2,rG);
        c = -.0001;
        while (delta < c && numEdgesCrossing>=1 && numEdgesCrossing_fiedler>=1)
            [posC,negC] =genCut(e2);
            sizeFiedCut = min(size(posC,2),size(negC,2));
            amountRemovedAcross = 0;
            [rG,amountRemovedAcross,amountRemovedTotal,amountAdded] = fiedRoundEigScore_2(rG,[delta],v,v,sum(rG),[budget],amountRemovedAcross,threshold,amountRemovedTotal);
            normLapY = computeNormalizedLaplacian(rG);
            [e2,lambda2new] = genSecondCutAndEig(normLapY);
            if abs(lambda2-lambda2new) < .000000000000001
                disp('did not make progress');
                break
            end
            [v,sizeCut,k] = findValidCut(.5*trueCutSize,normLapY);
            delta = targetSpec(k)-kthEig(normLapY,k);
            threshold = .5*abs(delta);
            numEdgesCrossing = compCrossEdges(v,rG);
            [e2,lambda2] = genSecondCutAndEig(normLapY);
            numEdgesCrossing_fiedler = compCrossEdges(e2,rG);
            count = count+1;
            if delta < 0
                budget = numEdgesCrossing*budgetParam;
            end
            clear normLapY
        end
    end



    

    function [F, G,e,tracePenalty,Y] = binary_2(X,gradFlag)
        Y = findAdj_linear(X);
        dY_ = sum(Y);
        dY = inv(diag(sqrt(abs(dY_))));
        entryPenalty = (Y.^2-2*Y.^3+Y.^4);
        F = 0;
        G=zeros(n);
        %Derivative of entryPenalty wrt to Y
        tempEntryPen = (4*Y.^3-6*Y.^2+2*Y);
        entryPenaltyDerivative = tempEntryPen-diag(diag(tempEntryPen));
        e = sum(sum(entryPenalty))-sum(diag(entryPenalty));
        F = F + e;
        tracePenalty = sum(diag(Y).^2);
        F = F + tracePenalty;        
        %Compute gradient if gradFlag is set
        if (gradFlag==1)
            G=compEntryPenaltyDerivative_tensor(X,entryPenaltyDerivative,zeros(n),Y);
            G = compTraceDerivative_fast(X,diag(diag(Y)),diag(sum(Y)),G);
        end
    end



    function [G] = compEntryPenaltyDerivative_tensor(X,entryDer,G,Y)
        dS = sqrt(diag(sum(Y)));
        G = G-2*(dS*entryDer)*(dS*X*Lambda);
    end

    function [G] = compTraceDerivative_fast(X,Ydiag,degree,G)
        G = G - 4*Ydiag*degree*X*Lambda;
    end


    function [NormLap] = computeNormalizedLaplacian(A)
        degree = sum(A);
        degree(degree==0)=1;
        degreeSeq = diag(degree);
        invSqrtDeg = diag(sqrt(diag(degreeSeq)).^-1);
        laplacian = degreeSeq - A;
        NormLap = invSqrtDeg*laplacian*invSqrtDeg;
    end


    function [cutGraph] = genFractional(Y_c)
        [cutGraph] = genFractionalLp(Y_c,sum(Y_c));
    end






    function [startX,spectrumConfigModel,configGraph,raw] = configurationModelBenchmark(degreeSequence)
        noZeroDeg = 0;
        while noZeroDeg == 0
            [configGraph,raw] = configurationModel(degreeSequence,floor(sum(degreeSequence)/2));
            if size(find(~(sum(configGraph))),2)==0
                noZeroDeg = 1;
            end
        end
        RandLap = computeNormalizedLaplacian(configGraph);
        [startX,Dthrowout] = eig(RandLap);
        spectrumConfigModel = diag(Dthrowout);
        spectrumConfigModel = sort(spectrumConfigModel);
    end





    


    function [A]=driveDownSpecGap(A,Spec)
        rng shuffle;
        zeroDegNodes = find(~(sum(A)));
        A(zeroDegNodes, :) = [];
        A(:,zeroDegNodes) = [];
        normLapY = computeNormalizedLaplacian(A);
        [Xy,prevspec,perm_1] = compSortedSpectrum(normLapY);
        [fiedCut_Vertices,fiedCutComp_Vertices] = genkthCut(Xy,perm_1,2);
        [e2,lambda2] = genSecondCutAndEig(normLapY);
        numCrossBip = compCrossEdges(e2,A);
        %size(A,1)
        A(A<.0000000000001)=0;
        not_succ=1;
        %largest connected component
        if(numCrossBip==0)
            B = largestcomponent(A);
            A = A(B,B);
            normLapAfterCutFix = computeNormalizedLaplacian(A);
            [Xy,prevspecNormLap,perm_1] = compSortedSpectrum(normLapAfterCutFix);
            [fiedCut_Vertices,fiedCutComp_Vertices] = genkthCut(Xy,perm_1,2);
            crossEdges = A(fiedCut_Vertices,fiedCutComp_Vertices);
            numCrossBip = sum(sum(crossEdges));
        end
        %if the fiedler cut is still an extrmeley unbalanced one, cut
        %the small part off and proceed
        while min(size(fiedCut_Vertices,2),size(fiedCutComp_Vertices,2)) < 5
        %while min(size(fiedCut_Vertices,2),size(fiedCutComp_Vertices,2)) < 5 & while min(size(fiedCut_Vertices,2),size(fiedCutComp_Vertices,2)) > 0
            %disp('small');
            if size(fiedCut_Vertices,2) < size(fiedCutComp_Vertices,2)
                A= A(fiedCutComp_Vertices,fiedCutComp_Vertices);
            else
                A= A(fiedCut_Vertices,fiedCut_Vertices);
            end
            normLapAfterCutFix = computeNormalizedLaplacian(A);
            [Xy,prevspecNormLap,perm_1] = compSortedSpectrum(normLapAfterCutFix);
            [fiedCut_Vertices,fiedCutComp_Vertices] = genkthCut(Xy,perm_1,2);
        end
        normLapY = computeNormalizedLaplacian(A);
        [Xy,prevspec,perm_1] = compSortedSpectrum(normLapY);
        [fiedCut_Vertices,fiedCutComp_Vertices] = genkthCut(Xy,perm_1,2);
        crossEdges = A(fiedCut_Vertices,fiedCutComp_Vertices);
        numCrossBip = sum(sum(crossEdges));
        [e2,lambda2] = genSecondCutAndEig(normLapY);
        while lambda2-targetSpec(2) > .0001 && numCrossBip >= 1 && not_succ == 1
            %disp('iteration of drive down spec gap');
            [A] = criticalCutRounding(A,Spec);
            A(A<.0000000000001)=0;
            normLap = computeNormalizedLaplacian(A);
            lambda2_prev = lambda2;
            [e2,lambda2] = genSecondCutAndEig(normLap);
            %disp('abs difference')
            %abs(lambda2-lambda2_prev)
            if abs(lambda2-lambda2_prev)<.000000001
                not_succ = 0;
            end
            if(numCrossBip==0)
                B = largestcomponent(A);
                A = A(B,B);
                normLapAfterCutFix = computeNormalizedLaplacian(A);
                [Xy,prevspecNormLap,perm_1] = compSortedSpectrum(normLapAfterCutFix);
                [fiedCut_Vertices,fiedCutComp_Vertices] = genkthCut(Xy,perm_1,2);
                crossEdges = A(fiedCut_Vertices,fiedCutComp_Vertices);
                numCrossBip = sum(sum(crossEdges));
            end
            %if the fiedler cut is still an extrmeley unbalanced one, cut
            %the small part off and proceed
            %while min(size(fiedCut_Vertices,2),size(fiedCutComp_Vertices,2)) < 5 & while min(size(fiedCut_Vertices,2),size(fiedCutComp_Vertices,2)) > 0
            while min(size(fiedCut_Vertices,2),size(fiedCutComp_Vertices,2)) < 5
                if size(fiedCut_Vertices,2) < size(fiedCutComp_Vertices,2)
                    A= A(fiedCutComp_Vertices,fiedCutComp_Vertices);
                else
                    A= A(fiedCut_Vertices,fiedCut_Vertices);
                end
                normLapAfterCutFix = computeNormalizedLaplacian(A);
                [Xy,prevspecNormLap,perm_1] = compSortedSpectrum(normLapAfterCutFix);
                [fiedCut_Vertices,fiedCutComp_Vertices] = genkthCut(Xy,perm_1,2);
                not_succ=1
            end
            normLap = computeNormalizedLaplacian(A);
            [e2,lambda2] = genSecondCutAndEig(normLap);
            [Xy,prevspecNormLap,perm_1] = compSortedSpectrum(normLap);
            [fiedCut_Vertices,fiedCutComp_Vertices] = genkthCut(Xy,perm_1,2);
            crossEdges = A(fiedCut_Vertices,fiedCutComp_Vertices);
            numCrossBip = sum(sum(crossEdges));
        end
        normLap = computeNormalizedLaplacian(A);
        [Xy,prevspecNormLap,perm_1] = compSortedSpectrum(normLap);
        [fiedCut_Vertices,fiedCutComp_Vertices] = genkthCut(Xy,perm_1,2);
        numCrossBip = sum(sum(A(fiedCut_Vertices,fiedCutComp_Vertices)));
        %round maximum value to 1 if number of edges crossing is smaller than 1
        while numCrossBip < 1
            m = max(max(A(fiedCut_Vertices,fiedCutComp_Vertices)));
            [locx,locy]=find(A==m);
            u = locx(1);
            v = locy(1);
            A(u,v)=1.0;
            A(v,u)=1.0;
            [fiedCut_Vertices,fiedCutComp_Vertices] = genkthCut(Xy,perm_1,2);
            crossEdges = A(fiedCut_Vertices,fiedCutComp_Vertices);
            numCrossBip = sum(sum(crossEdges));
        end
    end
%%%%%%%%%%Stiefel

    function [t] = lineSearchStep(X0, obj, objstep, dirDerivative, descentDirection, maxIter, stepMax, c1, c2)
        %X0 start point
        %obj is objective function we are trying to optimize over the
        %manifold
        %objstep evaluates the objective function by taking a step in the
        %descent direction
        %dirDerivative comptues the directional derivative at points along
        %the descent direction
        %maxIter is maximum number of iterations
        %stepMax is maximum step-size
        %will find a step size t \in (t0,tMax)
        t0 = 0;
        tMax = stepMax;
        %disp('evaluating no step');
        [F0,G0] = feval(obj,X0,1);
        nu = feval(descentDirection,X0,G0);
        FPrime0 = feval(dirDerivative,X0,t0,nu,G0);
        FPrev = F0;
        t = (t0+tMax)/2;
        tPrev = t0;
        for i = 1:maxIter
            [F,G] = feval(objstep,t,X0,nu,obj);
            if (F > F0 + c1*t*FPrime0 || (i >1 && F >= FPrev))
                t = zoom(tPrev,t,X0, F0, FPrime0, nu, objstep, dirDerivative,obj,c1,c2);
                break;
            end
            FPrimet = feval(dirDerivative,X0,t,nu,G);
            if (abs(FPrimet) <= -c2*FPrime0)
                break;
            end
            if (FPrimet >= 0 )
                t = zoom(t,tPrev, X0, F0, FPrime0, nu, objstep, dirDerivative,obj,c1,c2);
                break;
                
            end
            tPrev = t;
            FPrev = F;
            t = (t + tMax)/2;
        end
        if(i == maxIter)
            disp('maxing out');
        end
    end

    function [t] = zoom(tLo,tHi,X0,F0, FPrime0, nu, objstep, dirDerivative,obj,c1,c2)
        Flo = feval(objstep,tLo,X0,nu,obj);
        cap = 50;
        for i = 1:cap
            t = (tLo + tHi)/2;
            [F,G] = feval(objstep,t,X0,nu,obj);            
            if (F > F0 + c1*t*FPrime0 || F >= Flo)
                tHi = t;
            else
                FPrimet = feval(dirDerivative,X0,t,nu,G);
                if (abs(FPrimet) <= -c2*FPrime0)
                    break;
                end
                if (FPrimet *(tHi-tLo) >= 0)
                    tHi = tLo;
                end
                tLo = t;
                Flo = F;
            end
            if (abs(tHi-tLo) < .0000001)
                t=tLo;
                break;
            end
            
        end
        if (i == cap)
            disp('zoomed for 50');
            t = tLo;
        end
    end

     function [XP] = polarStep(X,nu,t)
        %EP is a point on the tangent plane
        EP = nu*t;
        %Computes point XP on the manifold by using the polar
        %decomposition retraction method, computed via
        %eigenvalue decomposition
        %disp('is real EP');
        %disp(isreal(EP));
        CP = eye(n) + EP'*EP;
        %disp(isreal(CP));
        [VP,DP] = eig(CP);
        %disp('taking inverse the problem');
        %disp(isreal(VP));
        %disp(isreal(inv(VP)));
        sCP = VP*sqrt(DP)*inv(VP);
        XP = (X+EP)*inv(sCP);
    end

    function [nu] = steepestDescentPolar(X,G)
        %computes gradient related steepest descent direction at X.
        %this computes a direction in the tangent space that moves in the
        %gradient of F restricted to the manifold
        nu = X*G'*X - G;
        %disp('descent direction real');
        %disp(isreal(nu));
    end

    function [F, G,e,tracePenalty,Y] = binary(X,gradFlag)
        %objective is to find X which makes template generated using
        %template gen LP the most binary
        %call findAdj to find template Y that is consistent with laplacian
        %XDX'.
        Y = genTemplate(X);
        dY_ = sum(Y);
        dY = inv(diag(sqrt(abs(dY_))));
        entryPenalty = (Y.^2-2*Y.^3+Y.^4);
        F = 0;
        G=zeros(n);
        %Derivative of entryPenalty wrt to Y
        tempEntryPen = (4*Y.^3-6*Y.^2+2*Y);
        entryPenaltyDerivative = tempEntryPen-diag(diag(tempEntryPen));
        %entry penalty
        e = sum(sum(entryPenalty))-sum(diag(entryPenalty));
        F = F + e;
        %add additional penalty to keep trace equal to zero
        tracePenalty = sum(diag(Y).^2);
        F = F + tracePenalty;
        %Compute gradient if gradFlag is set
        if (gradFlag==1)
            G=compEntryPenaltyDerivative_tensor(X,entryPenaltyDerivative,zeros(n),Y);
            G = compTraceDerivative_fast(X,diag(diag(Y)),diag(sum(Y)),G);
        end
    end


    function [F,G]= polarObj(t,X,nu,obj)
        %t is the step-size, X is where we are at, nu is the descent
        %direction, obj is the function I am trying to optimize
        %first compute the point on the manifold stepping along direction
        %nu in the tangent space
        XP = polarStep(X,nu,t);
        %now i can compute the objective value using the objective function
        [F,G]=feval(obj,XP,1);
    end

    function [dirDer] = polarDer(X,t,nu,G)
        %X is the point we are at
        %nu is the descent direction
        %t is the step-size
        %G is the gradient of F at X
        %this computes the direcitonal derivative of F at X+(nu)t in direction
        %nu with respect to t where. In order to do so, we have to compute
        %the gradient of F restricted to the manifold.
        XP = polarStep(X,nu,t);
        %disp('is step still real');
        %disp(isreal(XP))
        dirDer = trace(nu'*((eye(n)-XP*XP')*G + X*(.5*(XP'*G-G'*XP))));
    end

    function [XP] = cayleyStep(X,nu,t)
        %return a new point
        Z = inv(eye(n) + (t/2)*nu);
        %Y(tau)=XP is a trial orthogonal basis
        XP = Z*(eye(n) - (t/2)*nu)*X;
        
    end

    function [nu] = steepestDescentCayley(X,G)
        %return the cayley transform
        nu = G*X'-X*G';
    end

    function [F,G]= cayleyObj(t,X,nu,obj)
        XP = cayleyStep(X,nu,t);
        [F,G] = feval(obj, XP,1);
    end

    function [dirDer] = cayleyDer(X,t,nu,G)
        Z = inv(eye(n) + (t/2)*nu);
        XP = cayleyStep(X,nu,t);
        XPrime = -Z*nu*.5*(X+XP);
        %Derivative of F with respect to tau at Y(tau)=XP
        dirDer = trace(G'*XPrime);
    end
%
    function [X,newGraph] = lineSearch(startX,obj,valueStep, step, directionalDer, descent, maxIter, stepMax, c1,c2,numSteps)
        [F,G,e,t,Y] = obj(startX,1);
        for i = 1:numSteps
            [tstar]= lineSearchStep(startX,obj,valueStep, directionalDer, descent, maxIter, stepMax, c1,c2);
            nu = descent(startX,G);
            X = step(startX,nu,tstar);
            [FX,GX,entryPen,tracePen,curGraph] = obj(X,1);
            if (F - FX <= .000001)
                break;
            else
                startX=X;
                F=FX;
                G = GX;
            end
        end
        newGraph=curGraph;
    end





%%%%%%%%%%%%%%Gen Models%%%%%%%%%%%%
function [A,BeforeRound] = configurationModel(degreeSequence,m)
        rng shuffle
        A = zeros(n);
        stubs = [];
        for vertex = 1:n
            d = degreeSequence(vertex);
            for miniVertex = 1:d
                stubs = [stubs vertex];
            end
        end
        for i = 1:m
            s = size(stubs,2);
            uIndex = randi([1 s],1,1);
            u = stubs(uIndex);
            stubs(uIndex) = [];
            vIndex = randi([1 s-1], 1,1);
            v = stubs(vIndex);
            stubs(vIndex) = [];
            A(u,v) = A(u,v) + 1;
            A(v,u) = A(v,u) + 1;
        end
        BeforeRound = A;
        %disp('num multi');
        multi = 0;
        for i = 1:n
            for j = 1:n
                if (A(i,j) >1)
                    multi = multi+A(i,j)-1;
                    %re-wire
                    A(i,j)=1;
                    A(j,i)=1;
                end
            end
        end
        % disp(multi);
        A = A-diag(diag(A));
    end



%%%%%%%%%%%%%%%%%UTILS%%%%%%%%%%%%%%%%%%%%%
    function recordSpec(rG,graphFileName,specFileName,loc,numGraph)
        rG = rG - diag(diag(rG));
        zeroDegs=findZero(rG);
        %if a degree is zero, computation of symmetric normalized laplacian
        %will have a division by zero. Removes nodes from graph before
        %compuation
        if size(zeroDegs,2)>0
            zeroDegs = sort(zeroDegs);
            numZeroDegs = size(zeroDegs,2);
            for i = 1:numZeroDegs
                rG(zeroDegs(numZeroDegs+1-i),:)=[];
                rG(:,zeroDegs(numZeroDegs+1-i))=[];
            end
        end
        rG(isnan(rG)) = 0;
        rG(isinf(rG)) = 0;
        normLap = computeNormalizedLaplacian(rG);
        [x0, specNormLap] = eig(normLap);
        specNormLap = sort(diag(specNormLap));
        graphFileName = strcat(strcat(strcat(graphFileName,'_'),string(numGraph)),'.csv');
        specFileName = strcat(strcat(strcat(specFileName,'_'),string(numGraph)),'.txt');
        csvwrite(strcat(loc,graphFileName),rG);
        dlmwrite(strcat(loc,specFileName),real(specNormLap));
        clear normLap
    end
    


    function [A] = readGraph(head)
        graphLoc = strcat('graphs_gml/',head);
        graphFile = strcat(graphLoc,'.gml');
        %read in the graph
        graph = read_gml(graphFile);
        A = get_adjacency_matrix(graph);
        A = A + A';
        A = A - diag(diag(A));
        A(A>.01)=1;
        A(A<.01)=0;
    end

function [posC,negC] = genFiedCut(Xy,perm)
        Xy=Xy(:,perm);
        ey2 = Xy(:,2);
        posC = [];
        negC = [];
        for i = 1:n
            if ey2(i) > 0
                posC = [posC,i];
            else
                negC = [negC,i];
            end
        end
    end

    function [posC,negC] = genCut(ey2)
        posC = [];
        negC = [];
        for i = 1:size(ey2)
            if ey2(i) > 0
                posC = [posC,i];
            else
                negC = [negC,i];
            end
        end
    end

    function [numEdges] = compCrossEdges(e2,A)
        [posC,negC] = genCut(e2);
        crossEdges = A(posC,negC);
        numEdges = sum(sum(crossEdges))/2;
    end

    function [posC,negC] = genkthCut(Xy,perm,k)
        Xy=Xy(:,perm);
        ey2 = Xy(:,k);
        posC = [];
        negC = [];
        num = size(ey2,1);
        for i = 1:num
            if ey2(i) > 0
                posC = [posC,i];
            else
                negC = [negC,i];
            end
        end
    end

    function [Xy,sortedSpec,perm_1] = compSortedSpectrum(X)
        [Xy,curSpec] = eig(X);
        curSpec = diag(curSpec);
        [useless,perm_1]=sort(curSpec);
        sortedSpec=curSpec(perm_1);
    end

    function [ek] = kthEig(X,k)
        [Xy,curSpec] = eig(X);
        curSpec = diag(curSpec);
        [useless,perm_1]=sort(curSpec);
        sortedSpec=curSpec(perm_1);
        ek = sortedSpec(k);
    end

    

    function [e,lambda] = genSecondCutAndEig(X)
        [V,w] = eigs(X,2,'SR');
        if w(1,1)<w(2,2)
            lambda = w(2,2);
            e = V(:,2);
        else
            lambda = w(1,1);
            e = V(:,1);
        end
    end

     function [nonBinPairs] = compNonBinPairs(Y)
        nonBinPairs = [];
        numRows = size(Y,1);
        numCols = size(Y,2);
        for i = 1:numRows
            for j = i:numCols
                if (Y(i,j)<1 && Y(i,j) > 0.000000001)
                    nonBinPairs = [nonBinPairs; i,j];
                end
            end
        end
     end

    function [indices] = findZero(A)
        An = A - diag(diag(A));
        An(An>=1)=1;
        indices=find(sum(An)==0);
    end
end

