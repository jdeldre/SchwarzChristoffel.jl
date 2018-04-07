module Reindex


export ShiftReindex,OddReindex,Reflect

struct ShiftReindex{T,N,O}
  array :: Vector{T}
end

function ShiftReindex(v::Vector{T},o::Int) where T
  N = length(v)+o
  ShiftReindex{T,N,o}(v)
end
(vec::ShiftReindex{T,N,O})(i::Int) where {T,N,O} =
                    (i > N || i < 1+O) ? T(0) : vec.array[i-O]
(vec::ShiftReindex{T,N,O})(ir::AbstractArray{M}) where {T,N,O,M} =
                    (vec::ShiftReindex{T,N,O}).(collect(ir))

struct OddReindex{T,N}
   array :: Vector{T}
end
function OddReindex(v::Vector{T}) where T
    N = length(v)-1
    OddReindex{T,N}(v)
end
(vec::OddReindex{T,N})(i::Int) where {T,N} =
        abs(i) > N ? T(0) : i < 0 ? conj(vec.array[1-i]) : vec.array[1+i]
(vec::OddReindex{T,N})(ir::AbstractArray{M}) where {T,N,M} =
                    (vec::OddReindex{T,N}).(collect(ir))
ReindexedArray = Union{ShiftReindex,OddReindex}

struct Reflect{T}
   array :: T
end
(vec::Reflect{<:ReindexedArray})(i) = vec.array(-i)
(vec::Reflect{<:AbstractArray})(i)  = vec.array[-i]



end
