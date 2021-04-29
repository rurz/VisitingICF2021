using VisitingICF2021
using Images

function imgload(path::String)
    return load(path)
end

function chsplit(img, c)
    chs = Images.channelview(img)
    if c == 1
        return chs_r = chs[1,:,:]
    elseif c == 2
        return chs_g = chs[2,:,:]
    elseif c == 3
        return chs_b = chs[3,:,:]
    end
end
