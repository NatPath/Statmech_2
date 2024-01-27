import Pkg; Pkg.add("Plots")

#usings
using Plots


# make sudoko

box_index_from_grid_indices(i,j,l)=Int8(1+closest_n_multiple(i-1,l)+closest_n_multiple(j-1,l)/l)

grid_indices_from_box_index(n,l)=(Int8((n-1)÷l+1),Int8((n-1)%l+1))

closest_n_multiple(x,n)=x-x%n

struct Sudoku
    l::Int8
    board::Matrix{Int8}
    rows_histograms::Matrix{Int8}
    cols_histograms::Matrix{Int8}
    boxes_histograms::Matrix{Int8}
    function Sudoku(l::Int8,init_str::String)
        board=zeros(Int8,l^2,l^2)
        rows_histograms=zeros(Int8,l^2,l^2)
        cols_histograms=zeros(Int8,l^2,l^2)
        boxes_histograms=zeros(Int8,l^2,l^2)
        if init_str=="rows"
            for i in range(1,l^2)
                for j in range(1,l^2)
                    board[i,j]=i
                    cols_histograms[i,j]=1
                end
                #init p'th box histogram
                p=i
                (i′,j′)=grid_indices_from_box_index(p,l)
                for k in range(1,l)
                    boxes_histograms[p,closest_n_multiple(i′-1,l)+k]=l               
                end
                rows_histograms[i,i]=l^2
            end
            
        end
        if init_str=="columns"
            for i in range(1,l^2)
                for j in range(1,l^2)
                    board[i,j]=j
                    rows_histograms[i,j]=1
                end
                #init p'th box histogram
                p=i
                (i′,j′)=grid_indices_from_box_index(p,l)
                for k in range(1,l)
                    boxes_histograms[p,closest_n_multiple(j′-1,l)+k]=l                    
                end
                cols_histograms[i,i]=l^2
            end
        end
        if init_str=="box"
            for i in range(1,l^2)
                for j in range(1,l^2)
                    board[i,j]=box_index_from_grid_indices(i,j,l)
                end
                for k in range(0,l-1)
                    rows_histograms[i,Int8(closest_n_multiple(i-1,l)+k+1)]=l
                    cols_histograms[i,Int8(closest_n_multiple(i-1,l)/l+1+k*l)]=l
                end
                boxes_histograms[i,i]=l^2
            end
        end
        new(l,board,rows_histograms,cols_histograms,boxes_histograms)
    end
end

function countPairsInRow(sudoku::Sudoku, i)
    sum=0
    for val in range(1,sudoku.l^2)
        sum+=binomial(sudoku.rows_histograms[i,val],2)
    end
    return sum
end

function countPairsInColumn(sudoku::Sudoku, j)
    sum=0
    for val in range(1,sudoku.l^2)
        sum+=binomial(sudoku.cols_histograms[j,val],2)
    end
    return sum    
end

function countPairsInBox(sudoku::Sudoku, n)
    sum=0
    for val in range(1,sudoku.l^2)
        sum+=binomial(sudoku.boxes_histograms[n,val],2)
    end
    return sum 
end

function countPairsTotal(sudoku::Sudoku)
    sum=0
    for i in range(1,sudoku.l^2)
        sum+=countPairsInRow(sudoku,i)
    end
    for j in range(1,sudoku.l^2)
        sum+=countPairsInColumn(sudoku,j)
    end
    for n in range(1,sudoku.l^2)
        sum+=countPairsInBox(sudoku,n)
    end
    
    return sum 
end

function energy_of_sudoku(sudoku::Sudoku, highest_energy)
    return countPairsTotal(sudoku)/highest_energy
end

function pick_random_pair_of_sites(l)
    random_site=[rand(0:l^2-1),rand(0:l^2-1)]
    random_dist=[0,0]
    while (random_dist[1]==0 && random_dist[2]==0)
        random_dist[1]=rand(0:l^2-1)
        random_dist[2]=rand(0:l^2-1)
    end
    random_site′=((random_site+random_dist).%(l^2)).+1
    random_site=random_site.+1
    @assert (random_site[1]!=random_site′[1] || random_site[2]!=random_site′[2]) "Error: picked the same site $random_site $random_site′ $random_dist"
    return (random_site,random_site′)
end

function switch_sites!(sudoku,i,j,i′,j′)
    val=sudoku.board[i,j]
    val′=sudoku.board[i′,j′]
    if val==val′
        return
    end
    n=box_index_from_grid_indices(i,j,sudoku.l)
    n′=box_index_from_grid_indices(i′,j′,sudoku.l)
    sudoku.board[i,j]=val′
    sudoku.board[i′,j′]=val
    
    #update histograms
    sudoku.cols_histograms[j,val]-=1
    sudoku.cols_histograms[j,val′]+=1
    sudoku.rows_histograms[i,val]-=1
    sudoku.rows_histograms[i,val′]+=1
    sudoku.boxes_histograms[n,val]-=1
    sudoku.boxes_histograms[n,val′]+=1
    
    sudoku.cols_histograms[j′,val]+=1
    sudoku.cols_histograms[j′,val′]-=1
    sudoku.rows_histograms[i′,val′]-=1
    sudoku.rows_histograms[i′,val]+=1
    sudoku.boxes_histograms[n′,val]+=1
    sudoku.boxes_histograms[n′,val′]-=1
end

function acceptance_rate(E,E′,T)
    if ( E<E′ && T==0 )
        return 0
    end
    if E<E′
        return exp(-(E-E′)/T)
    end
    return 1
end

function print_histograms(sudoku)
    println("cols one by one:")
    for j in range(1,sudoku.l^2)
        println("col $j")
        println(sudoku.cols_histograms[j,:])
    end
    println("cols histogram together:")
    println(sudoku.cols_histograms)
    println("rows one by one:")
    for i in range(1,sudoku.l^2)
        println("row $i")
        println(sudoku.rows_histograms[i,:])
    end
    println("rows histogram together:")
    println(sudoku.rows_histograms)
    println("boxes one by one:")
    for n in range(1,sudoku.l^2)
        println("box $n")
        println(sudoku.boxes_histograms[n,:])
    end
    println("boxes histogram together:")
    println(sudoku.boxes_histograms)
end

function find_ground_states_sudoku(l,T_scale=1,boosting=0)
    sudoku=Sudoku(Int8(l),"box")
    total_num_of_pairs=countPairsTotal(sudoku)
    println("Total_num_of_pairs=$total_num_of_pairs")

    steps=10^(l+1+boosting)
    if T_scale==0
       temperatures=[0] 
    else
        temperatures=10*T_scale:-T_scale:0
    end
    energy=energy_of_sudoku(sudoku,total_num_of_pairs)
    switches=Array{Tuple}(undef,0)
    energies=[energy]
    for T in temperatures
        for n in range(1,steps)
            (site,site′)=pick_random_pair_of_sites(l)
            sites=(site,site′)
            i=site[1]
            j=site[2]
            i′=site′[1]
            j′=site′[2]
            last_energy=last(energies)
            switched=(i,j,i′,j′)
            switch_sites!(sudoku,i,j,i′,j′)
            energy′=energy_of_sudoku(sudoku,total_num_of_pairs)
            randnum=rand()
            if randnum>=acceptance_rate(last(energies),energy′,T)
                switch_sites!(sudoku,i,j,i′,j′)
                push!(energies,last(energies))
            else
                push!(energies,energy′)
                push!(switches,(i,j,i′,j′))
            end
        end
    end

    total_num_of_pairs=countPairsTotal(sudoku)
    println("Total_num_of_pairs=$total_num_of_pairs")
    ax=heatmap(sudoku.board, aspect_ratio=:equal, color=:viridis, clim=(0,l^2))
    res_plot=plot(0:steps*length(temperatures),energies)
    return (ax, total_num_of_pairs,res_plot)
end

res2_T0=find_ground_states_sudoku(2,0)

res2_T0_2=res2_T0[2]
println("ground state reached = $res2_T0_2")
res2_T0[1]

res2_T0[3]

res2=find_ground_states_sudoku(2)

res22=res2[2]
println("ground state reached = $res22")
res2[1]

res2[3]

res3_T0=find_ground_states_sudoku(3,0)

res3_T0_2=res3_T0[2]
println("ground state reached = $res3_T0_2")
res3_T0[1]

res3_T0[3]

res3=find_ground_states_sudoku(3)

res32=res3[2]
println("ground state reached = $res32")
res3[1]

res3[3]

res4_T0=find_ground_states_sudoku(4,0,1)

res4_T0_2=res4_T0[2]
println("ground state reached = $res4_T0_2")
res4_T0[1]

res4_T0[3]

res4′=find_ground_states_sudoku(4,1,1)

res4′2=res4′[2]
println("ground state reached = $res4′2")
res4′[1]

res4′[3]
