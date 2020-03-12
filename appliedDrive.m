function [consts,vecs,mats] = appliedDrive(consts,vecs,mats)

    mats.zxt = (exp(-(consts.YDirPresent*(abs(mats.yy)-min(abs(mats.yy),0.005+0*mats.yy)).^2+(mats.xx-mats.tt*(mats.v_b(1,1,consts.inflowCenter)).*consts.v_wave).^2)/consts.spread.^2/2)).*(mats.tt<=consts.signal_time).*(mats.tt>=0);

    mats.Fxt = 0*ndgrid(vecs.spacex, vecs.spacey, vecs.spacez, [0 0])+consts.F_0;
    
end

