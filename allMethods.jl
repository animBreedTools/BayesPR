
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

function stJWAS(phenoDataInRef::DataFrame,phenoDataInVal::DataFrame,genoData_All::DataFrame,trait::Int,BayesX::String,piValue::Float64,nChain::Int,nBurnin::Int,nThin::Int,varR::Float64,varG::Float64)
    #the changes (use of copy, and changes in the keywords of function) is because JWAS requires string. This changes the original data file if I use old versions of functions 
    phenoData_G4 = copy(phenoDataInRef)
    phenoData_G5 = copy(phenoDataInVal)
    genoData_Combined = copy(genoData_All)
    
    genoData_Combined[:ID]  = "ind".*string.(genoData_Combined[:ID])
    phenoData_G4[:ID] = "ind".*string.(phenoData_G4[:ID])
    phenoData_G5[:ID] = "ind".*string.(phenoData_G5[:ID])
    
    gInd      = genoData_Combined[:ID]
    gpInd     = intersect(genoData_Combined[:ID],phenoData_G4[:ID])
    gNoPInd   = setdiff(gInd,phenoData_G4[:ID])

    #not IDs, rows!
    refRows = [find(i -> i == j, phenoData_G4[:ID])[] for j in gpInd]
    phenoRef = phenoData_G4[refRows,:];

    #not IDs, rows!
    refRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gpInd]
    genoRef = genoData_Combined[refRows,:]; #first one is ID
    
    writecsv("refGeno",convert(Array,genoRef))

    model_equations = "pheno$trait = intercept" ;
    model1 = build_model(model_equations,varR);
    add_markers(model1,"refGeno",varG,separator=',',header=false);

    out = runMCMC(model1,phenoRef,Pi=piValue,estimatePi=false,chain_length=nChain,burnin=nBurnin,methods=BayesX,output_samples_frequency=nThin,output_file="MCMC_samples_$BayesX$(Int(piValue)).txt");
    #MCMC_marker_effects_output_file was changes to output_file

    #not IDs, rows!
    # first 200 is sires in G3 and G4 gNoPInd[401:end]
    testRows = [find(i -> i == j, phenoData_G5[:ID])[] for j in gNoPInd[401:end]]
    phenoTest = phenoData_G5[testRows,:];
    #not IDs, rows!
    # first 200 is sires gNoPInd[201:end]
    testRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gNoPInd[401:end]]
    genoTest = genoData_Combined[testRows,2:end];

    #ebvBayes = convert(Array{Int64},genoTest)*out["Posterior mean of marker effects"]
    snpEff   = convert(Array,readtable("MCMC_samples_$BayesX$(Int(piValue)).txt_marker_effects_pheno$trait.txt",separator=',',header=false))
    println("size of SNPEFF: $(size(snpEff))")
    snpEff   = mean(snpEff,dims=1)'
    ebvBayes = convert(Array{Int64},genoTest)*snpEff
    
    println("TRT $trait r in Tst ", cor(ebvBayes,convert(Array,phenoTest[Symbol("u$trait")])))
    r_Bayes = cor(ebvBayes,convert(Array,phenoTest[Symbol("u$trait")]))

    varE_Bayes = out["Posterior mean of residual variance"]
    
    varSNP_Bayes = convert(Array,readtable("MCMC_samples_$BayesX$(Int(piValue)).txt_marker_effects_variances.txt",header=false,separator=','))
    println("size of SNPEFF: $(size(varSNP_Bayes))")
    varSNP_Bayes = vcat(mean(varSNP_Bayes,dims=1)...)
    removeMe = "MCMC_samples_$BayesX$(Int(piValue)).txt_marker_effects_variances.txt"
    println("removeMe $removeMe removed")
    rm(removeMe)
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
    p    = mean(X,dims=1)./2.0
    X  .-= ones(Float64,nInd)*2p

    y    = convert(Array,phenoRef[Symbol("pheno$trait")])
    y    .-= mean(y)

    λ    = varR/varSNP

    βhat = X' * inv(X*X' + eye(nInd)*λ) * y

    #not IDs, rows!
    # first 200 in G3 and G4 are sires gNoPInd[201:end]
    testRows = [find(i -> i == j, phenoData_G5[:ID])[] for j in gNoPInd[401:end]]
    phenoTest = phenoData_G5[testRows,:];
    #not IDs, rows!
    # first 200 is sires gNoPInd[201:end]
    testRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gNoPInd[401:end]]
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
    p    = mean(X,dims=1)./2.0
    X  .-= ones(Float64,nInd)*2p

    y    = convert(Array,phenoRef[Symbol("pheno$trait")])
    y    .-= mean(y)

    Λi = inv(Λ)

    βhat = Λi*X' * inv(X*Λi*X' + eye(nInd)) * y

    #not IDs, rows!
    # first 200 in G3 and G4 are sires gNoPInd[201:end]
    testRows = [find(i -> i == j, phenoData_G5[:ID])[] for j in gNoPInd[401:end]]
    phenoTest = phenoData_G5[testRows,:];
    #not IDs, rows!
    # first 200 is sires gNoPInd[201:end]
    testRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gNoPInd[401:end]]
    genoTest = genoData_Combined[testRows,2:end];
    
    uHat = convert(Array,genoTest)*βhat

    r_wSNPBLUP = cor(uHat,convert(Array,phenoTest[Symbol("u$trait")]))
    return r_wSNPBLUP
end

function PBLUP(phenoData_G4::DataFrame,phenoData_G5::DataFrame,genoData_Combined::DataFrame,popPedigree::Array,trait::Int,varR::Float64,varG::Float64)
    gInd      = genoData_Combined[:ID]
    gpInd     = intersect(genoData_Combined[:ID],phenoData_G4[:ID])
    gNoPInd   = setdiff(gInd,phenoData_G4[:ID])

    #not IDs, rows!
    # first 200 in G3 and G4 are sires gNoPInd[201:end]
    testRows = [find(i -> i == j, phenoData_G5[:ID])[] for j in gNoPInd[401:end]]
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

function prepDataSSBR(phenoData_G4::DataFrame,genoData_Combined::DataFrame,popPedigree::Array,trait::Int)
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

function mmeSSBR(phenoData_G5::DataFrame,trait::Int,varSNP,varG,varR,Z1,X,X1,W,W1,y,y1,Ai11,J,M,nTot,gInd,ngInd,gNoPInd)    

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

    sol = gauss_seidel(lhs, rhs)
    #sol=lhs\rhs

    aHat  = J*sol[2] + M*sol[3:(length(sol)-n1)]  #BV from genotypes (either true or imputed)
    aHat[1:n1,:] += sol[(length(sol)-n1+1):end]   #this adds epsilon to ng individuals
    ebv = [[ngInd;gInd] aHat]

    testRows = [find(i -> i == j, ebv[:,1])[] for j in gNoPInd[401:end]];
    ebvPred = ebv[testRows,2]
    println("number of gNoPInd: $(length(testRows))")    
    testPhenoRows = [find(i -> i == j, phenoData_G5[:ID])[] for j in gNoPInd[401:end]];
    ebvTrue = phenoData_G5[testPhenoRows,Symbol("u$trait")]
    
    noPnoGInd = setdiff(phenoData_G5[:ID],gNoPInd[401:end])
    testRows2 = [find(i -> i == j, ebv[:,1])[] for j in noPnoGInd];
    ebvPred2 = ebv[testRows2,2]
    println("number of noPnoGInd: $(length(testRows2))")    
    testPhenoRows2 = [find(i -> i == j, phenoData_G5[:ID])[] for j in noPnoGInd];
    ebvTrue2 = phenoData_G5[testPhenoRows2,Symbol("u$trait")]
    
    r_ssSNPBLUP = [cor(ebvPred,ebvTrue) cor(ebvPred2,ebvTrue2)]
    return r_ssSNPBLUP
end

function mtJWAS(phenoDataInRef::DataFrame,phenoDataInVal::DataFrame,genoData_All::DataFrame,nTrait::Int,BayesX::String,piValue,nChain::Int,nBurnin::Int,nThin::Int,varR,varG)
    #the changes (use of copy, and changes in the keywords of function) is because JWAS requires string. This changes the original data file if I use old versions of functions 
    phenoData_G4 = copy(phenoDataInRef)
    phenoData_G5 = copy(phenoDataInVal)
    genoData_Combined = copy(genoData_All)
    
    genoData_Combined[:ID]  = "ind".*string.(genoData_Combined[:ID])
    phenoData_G4[:ID] = "ind".*string.(phenoData_G4[:ID])
    phenoData_G5[:ID] = "ind".*string.(phenoData_G5[:ID])

    gInd      = genoData_Combined[:ID]
    gpInd     = intersect(genoData_Combined[:ID],phenoData_G4[:ID])
    gNoPInd   = setdiff(gInd,phenoData_G4[:ID])

    #not IDs, rows!
    refRows = [find(i -> i == j, phenoData_G4[:ID])[] for j in gpInd]
    phenoRef = phenoData_G4[refRows,:];

    #not IDs, rows!
    refRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gpInd]
    genoRef = genoData_Combined[refRows,:]; #first one is ID

    writecsv("refGeno",convert(Array,genoRef))

    model_equations = "pheno1 = intercept
                       pheno2 = intercept";
    model1 = build_model(model_equations,varR);
    add_markers(model1,"refGeno",varG,separator=',',header=false);

    out = runMCMC(model1,phenoRef,Pi=piValue,estimatePi=false,chain_length=nChain,burnin=nBurnin,methods=BayesX,output_samples_frequency=nThin);

    #not IDs, rows!
    # first 200 is sires in G3 and G4 gNoPInd[401:end]
    testRows = [find(i -> i == j, phenoData_G5[:ID])[] for j in gNoPInd[401:end]]
    phenoTest = phenoData_G5[testRows,:];
    #not IDs, rows!
    # first 200 is sires gNoPInd[201:end]
    testRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gNoPInd[401:end]]
    genoTest = genoData_Combined[testRows,2:end];

    #ebvBayes = convert(Array{Int64},genoTest)*hcat(out["Posterior mean of marker effects"]...)
    snpEff   = [mean(convert(Array,readtable("MCMC_samples_marker_effects_pheno1.txt",separator=',',header=false))[Int(nBurnin/nThin)+1:end,:],dims=1)' mean(convert(Array,readtable("MCMC_samples_marker_effects_pheno2.txt",separator=',',header=false))[Int(nBurnin/nThin)+1:end,:],dims=1)']
    ebvBayes = convert(Array{Int64},genoTest)*snpEff

    println("r in Tst ", diag(cor(ebvBayes,convert(Array,phenoTest[[:u1,:u2]]))))
    r_Bayes =  diag(cor(ebvBayes,convert(Array,phenoTest[[:u1,:u2]])))

    varUhat = cov(ebvBayes)
    varE_Bayes = out["Posterior mean of residual variance"]
    
#    coVarSNP_Bayes = Array{Any}(0, 4)
    
    if BayesX=="BayesB"
        varData = CSV.read("MCMC_samples_marker_effects_variances.txt",delim=',',header=false)
        var1    = varData[collect(1:2:size(varData,1)),1]
        var2    = varData[collect(2:2:size(varData,1)),2]
        coVar12 = varData[collect(2:2:size(varData,1)),1]
        println("size of var1 $(size(var1))")
        var1    = reshape(var1,size(genoTest,2),Int(nChain/nThin))[:,Int(nBurnin/nThin)+1:end]   ##
        var2    = reshape(var2,size(genoTest,2),Int(nChain/nThin))[:,Int(nBurnin/nThin)+1:end]   ## Hao's JWAS prints out everything.
        coVar12 = reshape(coVar12,size(genoTest,2),Int(nChain/nThin))[:,Int(nBurnin/nThin)+1:end]## 
        println("size of var1 $(size(var1))")
        meanVar1    = mean(var1,dims=2)
        meanVar2    = mean(var2,dims=2)
        meanCoVar12 = mean(coVar12,dims=2)
#        coVarSNP_Bayes = vcat(coVarSNP_Bayes,[meanVar1 meanCoVar12 meanCoVar12 meanVar2])
        coVarSNP_Bayes = [meanVar1 meanCoVar12 meanCoVar12 meanVar2]
    elseif BayesX=="BayesC"
        varData = convert(Array,CSV.read("MCMC_samples_marker_effects_variances.txt",delim=',',header=false))[Int(nBurnin/nThin)+1:end,:]   ## Hao's JWAS prints out everything.
#        coVarSNP_Bayes = vcat(coVarSNP_Bayes,mean(varData,1))
        coVarSNP_Bayes = mean(varData,dims=1)
    end
#    removeMe = "MCMC_samples_$BayesX$(Int(piValue)).txt_variance.txt"
#    println("removeMe $removeMe removed")
#    rm(removeMe)
    return r_Bayes, varUhat, varE_Bayes, coVarSNP_Bayes
end

function prepDataSSBR_mt(phenoData_G4::DataFrame,genoData_Combined::DataFrame,popPedigree::Array,nTraits::Int)
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

    y_2Trait = Array{Float64}(undef,0,1)
    y1_2Trait = Array{Float64}(undef,0,1)
    Z11Z11 = Array{Float64}(undef,0,0)
    Z2Z2 = Array{Float64}(undef,0,0)
    XX   = Array{Float64}(undef,0,0)
    X1X1 = Array{Float64}(undef,0,0)
    WW   = Array{Float64}(undef,0,0)
    W1W1 = Array{Float64}(undef,0,0)
    JJ   = Array{Float64}(undef,0,0)

    for trait in 1:nTraits
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
    Z11 = Z1[:,1:n1]
            
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
   
    y_2Trait = vcat(y_2Trait,y)
    y1_2Trait = vcat(y1_2Trait,y1)
    Z11Z11 = cat([1,2],Z11Z11,Z11)
    Z2Z2 = cat([1,2],Z2Z2,Z2)
    XX   = cat([1,2],XX,X)
    X1X1 = cat([1,2],X1X1,X1)
    WW   = cat([1,2],WW,W)
    W1W1 = cat([1,2],W1W1,W1)
    JJ   = cat([1,2],JJ,J)
    end
    MM   = cat([1,2],M,M)
    return Z11Z11, XX, X1X1, WW, W1W1, y_2Trait, y1_2Trait, Ai11, JJ, MM, nTot, gInd, ngInd, gNoPInd 
end

function mmeSSBR_mt(phenoData_G5::DataFrame,nTraits::Int,coVarSNP,varG,varR,Z11Z11,XX,X1X1,WW,W1W1,y_2Trait,y1_2Trait,Ai11,JJ,MM,nTot,gInd,ngInd,gNoPInd)    

    n1 = length(ngInd)
    n2 = length(gInd)
    n3 = Int(size(MM,2)/nTraits)
    
    if size(coVarSNP,1)==1
        coVarSNP = vcat(fill(coVarSNP,n3)...)
        else
        coVarSNP = coVarSNP
    end        

    B = full([Diagonal(coVarSNP[:,1]) Diagonal(coVarSNP[:,2]);
    Diagonal(coVarSNP[:,3]) Diagonal(coVarSNP[:,4])])
    
    invR = inv(kron(varR,eye(Int(length(y_2Trait)/2)))) #assumes same number of pheno for each trait

    #invR1 = inv(kron(varR,eye(Int(length(y1_2Trait)/2)))) #assumes same number of pheno per trait
    
    lim1 = 1:Int(length(y1_2Trait)/2) #1:3000
    lim2 = Int(size(y_2Trait,1)/2)+1:Int(size(y_2Trait,1)/2)+Int(length(y1_2Trait)/2) #4001:7000
    invR1 = invR[vcat(lim1,lim2),vcat(lim1,lim2)] 
    
    C11 = XX'*invR*XX
    C12 = XX'*invR*WW
    C13 = X1X1'*invR1*Z11Z11
    C21 = WW'*invR*XX
    C22 = WW'*invR*WW + inv(B)
    C23 = W1W1'*invR1*Z11Z11
    C31 = Z11Z11'*invR1*X1X1
    C32 = Z11Z11'*invR1*W1W1
    C33 = Z11Z11'*invR1*Z11Z11+kron(inv(varG),Ai11);

    lhs = [C11 C12 C13;
           C21 C22 C23;
           C31 C32 C33];

    rhs = [XX'*invR*y_2Trait ; WW'*invR*y_2Trait ; Z11Z11'*invR1*y1_2Trait];

    sol = gauss_seidel(lhs, rhs)
    #sol = lhs\rhs
    
    aHat = JJ*sol[[2,4]] + MM*sol[5:(length(sol)-2*n1)]
    aHat = reshape(aHat,nTot,2)
    aHat[1:n1,:] += reshape(sol[(length(sol)-(2*n1)+1):end],size(ngInd,1),2)
    ebv = [[ngInd;gInd] aHat]
    
    testRows = [find(i -> i == j, ebv[:,1])[] for j in gNoPInd[401:end]];
    ebvPred = ebv[testRows,2:3]
    println("number of gNoPInd: $(length(testRows))")    
    testPhenoRows = [find(i -> i == j, phenoData_G5[:ID])[] for j in gNoPInd[401:end]];
    ebvTrue = phenoData_G5[testPhenoRows,[:u1,:u2]]

    noPnoGInd = setdiff(phenoData_G5[:ID],gNoPInd[401:end])
    testRows2 = [find(i -> i == j, ebv[:,1])[] for j in noPnoGInd];
    ebvPred2 = ebv[testRows2,2:3]
    println("number of noPnoGInd: $(length(testRows2))")    
    testPhenoRows2 = [find(i -> i == j, phenoData_G5[:ID])[] for j in noPnoGInd];
    ebvTrue2 = phenoData_G5[testPhenoRows2,[:u1,:u2]]
    
    r_ssSNPBLUP_mt = [diag(cor(ebvPred,convert(Array,ebvTrue))) diag(cor(ebvPred2,convert(Array,ebvTrue2)))]
    return r_ssSNPBLUP_mt
end

function runMTBayesPR(phenoDataInRef::DataFrame,phenoDataInVal::DataFrame,genoData_All::DataFrame,myMap,nChr::Int,rS::Int,nChain::Int,nBurnin::Int,nThin::Int,varR,varG)
    phenoData_G4 = copy(phenoDataInRef)
    phenoData_G5 = copy(phenoDataInVal)
    genoData_Combined = copy(genoData_All)
    
    genoData_Combined[:ID]  = "ind".*string.(genoData_Combined[:ID])
    phenoData_G4[:ID] = "ind".*string.(phenoData_G4[:ID])
    phenoData_G5[:ID] = "ind".*string.(phenoData_G5[:ID])

    gInd      = genoData_Combined[:ID]
    gpInd     = intersect(genoData_Combined[:ID],phenoData_G4[:ID])
    gNoPInd   = setdiff(gInd,phenoData_G4[:ID])

    #not IDs, rows!
    refRows = [find(i -> i == j, phenoData_G4[:ID])[] for j in gpInd]
    phenoRef = phenoData_G4[refRows,:];

    #not IDs, rows!
    refRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gpInd]
    genoRef = genoData_Combined[refRows,:]; #first one is ID
    
       #not IDs, rows!
    # first 200 is sires in G3 and G4 gNoPInd[401:end]
    testRows = [find(i -> i == j, phenoData_G5[:ID])[] for j in gNoPInd[401:end]]
    phenoTest = phenoData_G5[testRows,:];
    #not IDs, rows!
    # first 200 is sires gNoPInd[201:end]
    testRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gNoPInd[401:end]]
    genoTest = genoData_Combined[testRows,2:end];
    
    mtBayesPR(genoRef, phenoRef[[:pheno1, :pheno2]], myMap , nChr, rS, varG, varR, nChain, nBurnin, nThin, false)

    singleBeta = readtable("beta1Out"*"$rS",header=false)
    meanBeta1 = mean(convert(Array,singleBeta),1)'
    singleBeta = readtable("beta2Out"*"$rS",header=false)
    meanBeta2 = mean(convert(Array,singleBeta),1)'

    snpEff   = [meanBeta1 meanBeta2]
    ebvBayes = convert(Array{Int64},genoTest)*snpEff

    println("r in Tst ", diag(cor(ebvBayes,convert(Array,phenoTest[[:u1,:u2]]))))
    r_Bayes =  diag(cor(ebvBayes,convert(Array,phenoTest[[:u1,:u2]])))

    varUhat = cov(ebvBayes)

    covRegion = vcat(mean(convert(Array,readtable("covBetaOut"*"$rS",header=false)),dims=1)...)
    var1Region = vcat(mean(convert(Array,readtable("varBeta1Out"*"$rS",header=false)),dims=1)...)
    var2Region = vcat(mean(convert(Array,readtable("varBeta2Out"*"$rS",header=false)),dims=1)...)

    snpFile = readtable("snpInfo",header=false);
    regions = [searchsorted(snpFile[:x3],i) for i in 1:maximum(snpFile[:x3])];

    covSNP  = Array{Float64}(size(genoData[:,2:end],2))
    var1SNP = Array{Float64}(size(genoData[:,2:end],2))
    var2SNP = Array{Float64}(size(genoData[:,2:end],2))

    for i in 1:length(regions)
        covSNP[regions[i]]  .= covRegion[i] 
        var1SNP[regions[i]] .= var1Region[i] 
        var2SNP[regions[i]] .= var2Region[i] 
    end

    coVarSNP_BayesPR = [var1SNP covSNP covSNP var2SNP]

    p    = mean(convert(Array,genoRef[:,2:end]),1)./2.0


#    varG_BayesPR     = reshape(sum(coVarSNP_BayesPR.*(2*p.*(1-p))',1),2,2)
    varR_BayesPR = reshape(mean(convert(Array,readtable("varEOut"*"$rS",header=false)),1),2,2)
    gc()
    
    return r_Bayes, varUhat, varR_BayesPR, coVarSNP_BayesPR
end

function runSTBayesPR(phenoDataInRef::DataFrame,phenoDataInVal::DataFrame,genoData_All::DataFrame,trait,myMap,nChr::Int,rS::Int,nChain::Int,nBurnin::Int,nThin::Int,varR,varG)
    phenoData_G4 = copy(phenoDataInRef)
    phenoData_G5 = copy(phenoDataInVal)
    genoData_Combined = copy(genoData_All)
    
    genoData_Combined[:ID]  = "ind".*string.(genoData_Combined[:ID])
    phenoData_G4[:ID] = "ind".*string.(phenoData_G4[:ID])
    phenoData_G5[:ID] = "ind".*string.(phenoData_G5[:ID])

    gInd      = genoData_Combined[:ID]
    gpInd     = intersect(genoData_Combined[:ID],phenoData_G4[:ID])
    gNoPInd   = setdiff(gInd,phenoData_G4[:ID])

    #not IDs, rows!
    refRows = [find(i -> i == j, phenoData_G4[:ID])[] for j in gpInd]
    phenoRef = phenoData_G4[refRows,:];

    #not IDs, rows!
    refRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gpInd]
    genoRef = genoData_Combined[refRows,:]; #first one is ID
    
       #not IDs, rows!
    # first 200 is sires in G3 and G4 gNoPInd[401:end]
    testRows = [find(i -> i == j, phenoData_G5[:ID])[] for j in gNoPInd[401:end]]
    phenoTest = phenoData_G5[testRows,:];
    #not IDs, rows!
    # first 200 is sires gNoPInd[201:end]
    testRows = [find(i -> i == j, genoData_Combined[:ID])[] for j in gNoPInd[401:end]]
    genoTest = genoData_Combined[testRows,2:end];

    bayesPR(genoRef, phenoRef[Symbol("pheno$trait")], myMap , nChr, rS, varG, varR, nChain, nBurnin, nThin, false)

    meanBeta = readtable("betaOut"*"$rS",header=false)
    meanBeta   = mean(convert(Array,meanBeta),1)'

    ebvBayes = convert(Array{Int64},genoTest)*meanBeta
    
    println("r in Tst ", diag(cor(ebvBayes,convert(Array,phenoTest[Symbol("pheno$trait")]))))
    r_Bayes =  cor(ebvBayes,convert(Array,phenoTest[Symbol("pheno$trait")]))

    varUhat = cov(ebvBayes)

    varRegion = vcat(mean(convert(Array,readtable("varBetaOut"*"$rS",header=false)),dims=1)...)

    snpFile = readtable("snpInfo",header=false);
    
    regions = [searchsorted(snpFile[:x3],i) for i in 1:maximum(snpFile[:x3])];

    varSNP = Array{Float64}(size(genoData[:,2:end],2))

    for i in 1:length(regions)
        varSNP[regions[i]] .= varRegion[i] 
    end

    varR = mean(convert(Array,readtable("varEOut"*"$rS",header=false)),1)[]
    gc()
    
    return r_Bayes, varUhat, varR, varSNP
end


