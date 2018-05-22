
# adapted from http://morotalab.org/Mrode2005/relmat/createA.txt
function makeA(s::Array, d::Array)
    n = length(s)
    N = n + 1
    A = zeros(N, N)
    s = (s .== 0)*N + s
    d = (d .== 0)*N + d
for i in 1:n
    A[i,i] = 1.0 + A[s[i], d[i]]/2.0
        for j in (i+1):n
            if j > n break end
                A[i,j] = ( A[i, s[j]] + A[i,d[j]] )/2.0
                A[j,i] = A[i,j]
    end
    end
return(A[1:n, 1:n])
end

function stJWAS(phenoData_G4::DataFrame,phenoData_G5::DataFrame,genoData_Combined::DataFrame,trait::Int,BayesX::String,π::Float64,nChain::Int,nThin::Int,varR::Float64,varG::Float64)
    gInd      = genoData_Combined[:ID]
    gpInd     = intersect(genoData_Combined[:ID],phenoData_G4[:ID])
    gNoPInd   = setdiff(gInd,phenoData_G4[:ID])

    #not IDs, rows!
    refRows = [find(i -> i == j, phenoData_G4[:ID])[] for j in gpInd]
    phenoRef = phenoData_G4[refRows,:];

    #not IDs, rows!
    refRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gpInd]
    genoRef = genoData_Combined[refRows,2:end]; #first one is ID

    model_equations = "pheno$trait = intercept" ;
    model1 = build_model(model_equations,varR);
    add_markers(model1,genoRef,varG,separator=' ',header=true);

    out = runMCMC(model1,phenoRef,Pi=π,estimatePi=false,chain_length=nChain, methods=BayesX,output_samples_frequency=nThin,MCMC_marker_effects_file="MCMC_samples_$BayesX$(Int(π)).txt");

    #not IDs, rows!
    # first 200 is sires gNoPInd[201:end]
    testRows = [find(i -> i == j, phenoData_G5[:ID])[] for j in gNoPInd[201:end]]
    phenoTest = phenoData_G5[testRows,:];
    #not IDs, rows!
    # first 200 is sires gNoPInd[201:end]
    testRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gNoPInd[201:end]]
    genoTest = genoData_Combined[testRows,2:end];

    ebvBayes = convert(Array{Int64},genoTest)*out["Posterior mean of marker effects"]
    println("TRT $trait r in Tst ", cor(ebvBayes,convert(Array,phenoTest[Symbol("u$trait")])))
    r_Bayes = cor(ebvBayes,convert(Array,phenoTest[Symbol("u$trait")]))

    varE_Bayes = mean(out["MCMC samples for residual variance"])

    varSNP_Bayes = vcat(mean(convert(Array,readtable("MCMC_samples_$BayesX$(Int(π)).txt_variance.txt",header=false,separator=' ')),1)...)
    return r_Bayes, varE_Bayes, varSNP_Bayes
end

function SNPBLUP(phenoData_G4::DataFrame,phenoData_G5::DataFrame,genoData_Combined::DataFrame,trait::Int,varR::Float64,varSNP::Float64)
    gInd      = genoData_Combined[:ID]
    gpInd     = intersect(genoData_Combined[:ID],phenoData_G4[:ID])
    gNoPInd   = setdiff(gInd,phenoData_G4[:ID])

    #not IDs, rows!
    refRows = [find(i -> i == j, phenoData_G4[:ID])[] for j in gpInd]
    phenoRef = phenoData_G4[refRows,:];

    #not IDs, rows!
    refRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gpInd]
    genoRef = genoData_Combined[refRows,2:end]; #first one is ID

    nInd = size(genoRef,1)

    X    = convert(Array{Float64},genoRef)
    p    = mean(X,1)./2.0
    X  .-= ones(Float64,nInd)*2p

    y    = convert(Array,phenoRef[Symbol("pheno$trait")])
    y    .-= mean(y)

    λ    = varR/varSNP

    βhat = X' * inv(X*X' + eye(nInd)*λ) * y

    #not IDs, rows!
    # first 200 is sires gNoPInd[201:end]
    testRows = [find(i -> i == j, phenoData_G5[:ID])[] for j in gNoPInd[201:end]]
    phenoTest = phenoData_G5[testRows,:];
    #not IDs, rows!
    # first 200 is sires gNoPInd[201:end]
    testRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gNoPInd[201:end]]
    genoTest = genoData_Combined[testRows,2:end];

    uHat = convert(Array,genoTest)*βhat

    r_SNPBLUP = cor(uHat,convert(Array,phenoTest[Symbol("u$trait")]))
    return r_SNPBLUP
end

function wSNPBLUP(phenoData_G4::DataFrame,phenoData_G5::DataFrame,genoData_Combined::DataFrame,trait::Int,varR::Float64,varSNP::Array)
    gInd      = genoData_Combined[:ID]
    gpInd     = intersect(genoData_Combined[:ID],phenoData_G4[:ID])
    gNoPInd   = setdiff(gInd,phenoData_G4[:ID])

    #not IDs, rows!
    refRows = [find(i -> i == j, phenoData_G4[:ID])[] for j in gpInd]
    phenoRef = phenoData_G4[refRows,:];

    #not IDs, rows!
    refRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gpInd]
    genoRef = genoData_Combined[refRows,2:end]; #first one is ID

    varSNP =   full(Diagonal(varSNP))
    Λ      = full(Diagonal(varR./varSNP));

    nInd = size(genoRef,1)

    X    = convert(Array{Float64},genoRef)
    p    = mean(X,1)./2.0
    X  .-= ones(Float64,nInd)*2p

    y    = convert(Array,phenoRef[Symbol("pheno$trait")])
    y    .-= mean(y)

    Λi = inv(Λ)

    βhat = Λi*X' * inv(X*Λi*X' + eye(nInd)) * y

    #not IDs, rows!
    # first 200 is sires gNoPInd[201:end]
    testRows = [find(i -> i == j, phenoData_G5[:ID])[] for j in gNoPInd[201:end]]
    phenoTest = phenoData_G5[testRows,:];
    #not IDs, rows!
    # first 200 is sires gNoPInd[201:end]
    testRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gNoPInd[201:end]]
    genoTest = genoData_Combined[testRows,2:end];
    
    uHat = convert(Array,genoTest)*βhat

    r_wSNPBLUP = cor(uHat,convert(Array,phenoTest[Symbol("u$trait")]))
    return r_wSNPBLUP
end

function PBLUP(phenoData_G4::DataFrame,phenoData_G5::DataFrame,genoData_Combined::DataFrame,popPedigree::DataFrame,trait::Int,varR::Float64,varG::Float64)
    gInd      = genoData_Combined[:ID]
    gpInd     = intersect(genoData_Combined[:ID],phenoData_G4[:ID])
    gNoPInd   = setdiff(gInd,phenoData_G4[:ID])

    #not IDs, rows!
    # first 200 is sires gNoPInd[201:end]
    testRows = [find(i -> i == j, phenoData_G5[:ID])[] for j in gNoPInd[201:end]]
    phenoTest = phenoData_G5[testRows,:];

    nTot = size(popPedigree,1)
    allInd = collect(1:nTot)
    
    y   = fill(-9999.0,nTot)
    y[phenoData_G4[:ID,]] = phenoData_G4[Symbol("pheno$trait"),]
    
    y = y[find(y.!=-9999.0)]
    Z = eye(nTot)
    
    Z[:,setdiff(allInd,phenoData_G4[:ID])] .= 0
    Z = Z[find(sum(Z,2).!=0),:]


    
    G = (popPedigree*varG)
    R = varR*eye(length(y));
    
    println("sizeG $(size(G)) sizeR $(size(R))")

    V = Z*G*Z'+ R
    
    uHatPheno = G*Z'*inv(V)*(y.-mean(y))

    r_PBLUP = cor(uHatPheno[end-(size(phenoTest,1)-1):end],convert(Array,phenoTest[Symbol("u$trait")]))
    return r_PBLUP
end

function prepDataSSBR(phenoData_G4::DataFrame,genoData_Combined::DataFrame,popPedigree::DataFrame,trait::Int)
    nTot = size(popPedigree,1)
    allInd    = collect(1:nTot)
    gpInd     = intersect(genoData_Combined[:ID],phenoData_G4[:ID]);
    gInd      = genoData_Combined[:ID]
    ngInd     = setdiff(allInd,gInd)
    gNoPInd   = setdiff(gInd,phenoData_G4[:ID])
    pNoGInd   = setdiff(phenoData_G4[:ID],gInd)

    n1 = length(ngInd)
    n2 = length(gInd)
    
    popPedigree = popPedigree[[ngInd;gInd],[ngInd;gInd]]
    popPedigree[1:10,:]

    Ai = inv(popPedigree)

    Ai11 = Ai[1:n1,1:n1];
    Ai12 = Ai[1:n1,(n1+1):nTot];
    print(size(Ai11)," ", size(Ai12))

    n2, nCols = size(genoData_Combined)
    nMarkers  = nCols - 1
    M2 = convert(Array{Float32},genoData_Combined[:,2:end])
    print(size(M2))

    M1 = -Ai11\(Ai12*M2)
    M  = [M1;M2]

    yTemp   = fill(-9999.0,nTot)
    yTemp[phenoData_G4[:ID,]] = phenoData_G4[Symbol("pheno$trait"),]
    y = yTemp[[ngInd;gInd]]
    y1Temp = y[1:n1]
    y2Temp = y[(n1+1):end]
    y1 = y1Temp[y1Temp.!=-9999.0]
    y2 = y2Temp[y2Temp.!=-9999.0]
    y  = [y1;y2]

    Z2 = eye(nTot)
    Z2[:,gNoPInd] .= 0
    Z2 = Z2[gpInd,[ngInd;gInd]]

    Z1 = full(Diagonal((y1Temp.!=-9999.0)*1))
    Z1 = Z1[find(sum(Z1,2).!=0),:]
    Z1 = [Z1 zeros(length(pNoGInd), length(gInd))]

    n1 = size(M1,1)
    n2 = size(M2,1)
    J2 = -ones(n2,1)
    J1 = -Ai11\(Ai12*J2)
    J = [J1;J2]

    X  = [hcat(ones(n1), J1);
          hcat(ones(n2), J2)]
    X1 = Z1*X
    X2 = Z2*X
    X = [X1;X2]

    W1 = Z1*M
    W2 = Z2*M
    W  = [W1;W2]
    return Z1, X, X1, W, W1, y, y1, Ai11, J, M, nTot, gInd, ngInd, gNoPInd 
end

function mmeSSBR(phenoData_G5::DataFrame,trait::Int,varSNP,Z1,X,X1,W,W1,y,y1,Ai11,J,M,nTot,gInd,ngInd,gNoPInd)    

    n1 = length(ngInd)
    n2 = length(gInd)
    n3 = size(M,2)
    
    if length(varSNP)==1
        covarSNP = fill(varSNP,n3)
        else
        covarSNP = varSNP
    end
        
    D      = full(Diagonal(varR./covarSNP));

    λ1 = varR/varG

    Z11 = Z1[:,1:n1]
    lhs = [X'X     X'W             X1'Z11;
           W'X     W'W+D           W1'Z11;
           Z11'X1  Z11'W1          Z11'Z11+Ai11*λ1]
    rhs = [X'y; W'y; Z11'y1];

    sol=lhs\rhs

    aHat  = J*sol[2] + M*sol[3:(length(sol)-n1)]
    size(aHat)
    aHat[1:n1,:] += sol[(length(sol)-n1+1):end];
    ebv = [[ngInd;gInd] aHat]

    testRows = [find(i -> i == j, ebv[:,1])[] for j in gNoPInd];
    ebvPred = ebv[testRows,2]
    size(testRows)

    r_ssSNPBLUP = cor(ebvPred[201:end],phenoData_G5[Symbol("u$trait")])
    return r_ssSNPBLUP
end
