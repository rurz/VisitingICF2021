using VisitingICF2021
using Images

export imgload, chsplit

function imgload(path::String)
    return load(path)
end

function chsplit(img, c)
    chs = Images.channelview(img)
    if c == 1
        return chs[1,:,:]
    elseif c == 2
        return chs[2,:,:]
    elseif c == 3
        return chs[3,:,:]
    end
end
