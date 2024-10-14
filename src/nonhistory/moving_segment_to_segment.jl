@with_kw struct MovingSegmentToSegment{T <: Number} <: Setup @deftype T
    # Assuming movement in the x-direction
    x   
    y
    v
    D1
    H1
    D2
    H2
    σ = sqrt(x^2 + y^2)
end
