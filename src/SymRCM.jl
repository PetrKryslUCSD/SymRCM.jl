module SymRCM


# (C) 2020-2023, Petr Krysl


export symrcm


using Base: oneto
using SparseArrays
using SparseArrays: getcolptr


############################
# Abstract Graph Interface #
############################


const AbstractGraph{V} = Union{SparseMatrixCSC{<:Any, V}, AbstractVector{<:AbstractVector{V}}}


function neighbors(graph::SparseMatrixCSC, v::Integer)
    @view rowvals(graph)[nzrange(graph, v)]
end


function neighbors(graph::AbstractVector{<:AbstractVector}, v::Integer)
    graph[v]
end


function degree(graph::SparseMatrixCSC, v::Integer)
    getcolptr(graph)[v + 1] - getcolptr(graph)[v]
end


function degree(graph::AbstractVector{<:AbstractVector{V}}, v::Integer) where V
    n::V = length(graph[v])
    n
end


function Δ(graph::AbstractGraph{V}) where V
    maximum(vertices(graph); init=zero(V)) do v
        degree(graph, v)
    end
end


function nv(graph::SparseMatrixCSC{<:Any, V}) where V
    n::V = size(graph, 2)
    n
end


function nv(graph::AbstractVector{<:AbstractVector{V}}) where V
    n::V = length(graph)
    n
end


function vertices(graph::AbstractGraph)
    oneto(nv(graph))
end


function copygraph(graph::SparseMatrixCSC)
    copy(graph)
end


function copygraph(graph::AbstractVector{<:AbstractVector})
    deepcopy(graph)
end


###################################
# Reverse Cuthill-Mckee Algorithm #
###################################


"""
    symrcm(graph; sortbydeg::Bool=true) 

Reverse Cuthill-McKee node-renumbering algorithm.

- `sortbydeg`: Should the neighbor lists be sorted by column degree? The default is
  `true`, but often results of very similar quality are obtained when this is
  set to `false` and the lists are not sorted. The second option can be much
  faster, as the sorting is expensive when the neighbor lists are long.
"""
function symrcm(graph; sortbydeg::Bool=true)
    symrcm(graph, sortbydeg)
end


# Apply the reverse Cuthill-Mckee algorithm to each connected component of a graph.
function symrcm(graph::AbstractGraph{V}, sortbydeg::Bool) where V
    if sortbydeg
        graph = copygraph(graph)
        
        # sort neighbors
        scratch = Vector{V}(undef, Δ(graph))

        for j in vertices(graph)
            sort!(neighbors(graph, j); by=i -> degree(graph, i), scratch)
        end
    end

    symrcm_sorted(graph)
end


function symrcm_sorted(graph::AbstractGraph{V}) where V
    label = fill(false, nv(graph))
    order = Vector{V}(undef, nv(graph))

    if !iszero(nv(graph))
        # find psuedo-peripheral vertex
        root = argmin(vertices(graph)) do i
            degree(graph, i)
        end

        # compute Cuthill-Mckee ordering
        component, queue = bfs!(label, order, graph, root)

        for j in vertices(graph)
            if !label[j]
                # compute connected component
                component, queue = bfs!(label, queue, graph, j)

                # find psuedo-peripheral vertex
                root = argmin(component) do i
                    degree(graph, i)
                end

                # reset labels
                label[component] .= false

                # compute Cuthill-Mckee ordering
                bfs!(label, component, graph, root)
            end
        end
    end

    # reverse ordering    
    reverse!(order)
end


# Perform a breadth-first search of a graph.
@views function bfs!(label::AbstractVector{Bool}, queue::AbstractVector{V}, graph::AbstractGraph{V}, root::V) where V
    i = j = firstindex(queue)
    label[root] = true
    queue[j] = root

    @inbounds while i <= j
        for v in neighbors(graph, queue[i])
            if !label[v]
                j += 1
                label[v] = true
                queue[j] = v
            end
        end

        i += 1
    end

    queue[begin:j], queue[i:end]
end


############
# Old Code #
############


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


function symrcm(adjgr::Vector{Vector{T}}, degrees::Vector{T}) where T
    symrcm_sorted(adjgr)
end


end # module
