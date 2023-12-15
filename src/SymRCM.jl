module SymRCM

# (C) 2020-2023, Petr Krysl

export symrcm

using SparseArrays

"""
    adjgraph(A)

Compute the adjacency graph from a sparse matrix. 

The sparse matrix `A` is assumed to be symmetric.
The results will be wrong if it isn't.

- `sortbydeg`: Should the neighbor lists be sorted by column degree? The default is
  `true`, but often results of very similar quality are obtained when this is
  set to `false` and the lists are not sorted. The second option is much
  faster, as the sorting is expensive.
"""
function adjgraph(A::SparseMatrixCSC; sortbydeg = true)
    colptr = A.colptr
    rowval = A.rowval
    ncols = length(colptr)-1
    neighbors = Vector{Vector{eltype(colptr)}}(undef, ncols)
    cdeg = diff(colptr) # the degree is colptr[j+1]-colptr[j]
    for j in 1:ncols
        cstart = colptr[j]
        jdeg = cdeg[j]
        neighbors[j] = [rowval[cstart+m-1] for m in 1:jdeg]
    end
    # All of these sorts can be done in parallel,  they are totally independent.
    # The question is when to switch over to parallel execution so as to
    # amortize the cost of starting up threads.
    if sortbydeg
        for j in 1:ncols
            sort!(neighbors[j], by = j -> cdeg[j])
        end
    end
    return neighbors
end

"""
    adjgraph(conn, nfens)

Compute the adjacency graph from the array of connectivities of the elements
in the mesh.

# Examples
```
conn = [9 1 8 4;
       1 3 2 8;
       8 2 7 5;
       2 6 7 7];
nfens = 9;
adjgraph(conn, nfens)
```
should produce
```
9-element Array{Array{Int64,1},1}:
 [9, 8, 4, 3, 2]
 [1, 3, 8, 7, 5, 6]
 [1, 2, 8]
 [9, 1, 8]
 [8, 2, 7]
 [2, 7]
 [8, 2, 5, 6]
 [9, 1, 4, 3, 2, 7, 5]
 [1, 8, 4]
 ```
"""
function adjgraph(conn, nfens)
    neighbors = fill(Int[], nfens)
    for i = 1:length(neighbors)
        neighbors[i] = Int[]
    end
    for k in 1:size(conn, 1)
        for node1 in conn[k, :]
            for node2 in conn[k, :]
                if node1 != node2
                    push!(neighbors[node1], node2)
                end
            end
        end
    end
    # All of these can be done in parallel,  they are totally independent. The
    # question is when to switch over to parallel execution amortize the cost
    # of starting up threads.
    for i in 1:length(neighbors)
        neighbors[i] = unique(neighbors[i])
    end
    return neighbors
end

"""
    nodedegrees(adjgr::Vector{Vector{Int}})

Compute the degrees of the nodes in the adjacency graph.

conn = [9 1 8 4;
       1 3 2 8;
       8 2 7 5;
       2 6 7 7];
nfens = 9;
adjgr = adjgraph(conn, nfens)
nodedegrees(adjgr)

julia> degrees = node_degrees(adjgr)
9-element Array{Int64,1}:
 5
 6
 3
 3
 3
 2
 4
 7
 3
"""
function nodedegrees(adjgr::Vector{Vector{T}}) where {T}
    degrees = fill(zero(T), length(adjgr))
    for k in 1:length(adjgr)
        degrees[k] = length(adjgr[k])
    end
    return degrees
end

"""
    symrcm(adjgr::Vector{Vector{T}}, degrees::Vector{T}) where {T}

Reverse Cuthill-McKee node-renumbering algorithm.
"""
function symrcm(adjgr::Vector{Vector{T}}, degrees::Vector{T}) where {T}
    @assert length(adjgr) == length(degrees)
    # Initialization
    n = length(adjgr)
    ndegperm = sortperm(degrees) # sorted nodal degrees
    inR = fill(false, n) # Is a node in the result list?
    inQ = fill(false, n) # Is a node in the queue?
    R = T[]
    sizehint!(R, n)
    Q = T[] # Node queue
    sizehint!(Q, n)
    while true
        P = zero(T) # Find the next node to start from
        while !isempty(ndegperm)
            i = popfirst!(ndegperm)
            if !inR[i]
                P = i
                break
            end
        end
        if P == zero(T)
            break # That was the last node
        end
        # Now we have a node to start from: put it into the result list
        push!(R, P); inR[P] = true
        empty!(Q) # empty the queue
        append!(Q, adjgr[P]); inQ[adjgr[P]] .= true # put adjacent nodes in queue
        while length(Q) >= 1
            C = popfirst!(Q) # child to put into the result list
            inQ[C] = false # make note: it is not in the queue anymore
            if !inR[C]
                push!(R, C); inR[C] = true
            end
            for i in adjgr[C] # add all adjacent nodes into the queue
                if (!inR[i]) && (!inQ[i]) # contingent on not being in result/queue
                    push!(Q, i); inQ[i] = true
                end
            end
        end
    end
    return reverse(R) # reverse the result list
end

"""
    symrcm(A::SparseMatrixCSC; sortbydeg = true) 

Reverse Cuthill-McKee node-renumbering algorithm.

Compute the adjacency graph from a sparse matrix. The sparse matrix `A` is
assumed to be symmetric. The results will be wrong if it isn't.

- `sortbydeg`: Should the neighbor lists be sorted by column degree? The default is
  `true`, but often results of very similar quality are obtained when this is
  set to `false` and the lists are not sorted. The second option can be much
  faster, as the sorting is expensive when the neighbor lists are long.
"""
function symrcm(A::SparseMatrixCSC; sortbydeg = true) 
    ag = adjgraph(A; sortbydeg = sortbydeg)
    nd = nodedegrees(ag)
    return symrcm(ag, nd)
end

end # module
