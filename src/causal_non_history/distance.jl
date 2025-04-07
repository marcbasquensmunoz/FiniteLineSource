
compute_distance_3D(s, t) = sqrt((s.x - t.x)^2 + (s.y - t.y)^2 + (s.z - t.z)^2)
compute_distance_2D(s, t) = sqrt((s.x - t.x)^2 + (s.y - t.y)^2)

minimum_distance(source::PointSource, target::PointSource) = compute_distance_3D(source, target)

function minimum_distance(source::LineSource, target::LineSource)
    σ = compute_distance_2D(source, target)
    if target.D > source.D + source.H
        dist = sqrt(σ^2 + (target.D - source.D - source.H)^2)
    elseif source.D > target.D + target.H
        dist = sqrt(σ^2 + (source.D - target.D - target.H)^2)
    else 
        dist = σ
    end
    return dist
end
