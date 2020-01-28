% 	Written by Oliver T. Schmidt, 2012
%	University of Stuttgart
% 	E-mail: iagoschm@iag.uni-stuttgart.de
%
% 	2012/04 Oliver T. Schmidt           initial upload
%   2012/08 H. Kurz                     kennsatz-structure

function [ks, data] = eas2mat(file, ATTRLEN, UDEFLEN) 

fid=fopen(file);

ks.ATTRLEN=ATTRLEN;
ks.UDEFLEN=UDEFLEN;

%     identifyer: 20 byte character, e.g. "EAS3_I8R8 "
ks.identifyer=char(fread(fid, 20, 'char*1', 'ieee-be')');

%     file type: 8 byte integer (EAS2=1 or EAS3=2)
ks.dform=fread(fid, 1, 'int64', 'ieee-be')';

%     bform: 8 byte integer (IEEES=1, IEEED=2 or IEEEQ=3 for single double or quadruple bform, repsectively)
ks.bform=fread(fid, 1, 'int64', 'ieee-be')';

%     array sizes: each 8 byte integer
%         number of timesteps: nzs
ks.nzs=fread(fid, 1, 'int64', 'ieee-be')';

%         number of parameters: npar
ks.npar=fread(fid, 1, 'int64', 'ieee-be')';

%         number of gridpoints along dimension 1: ndim1
ks.ndim1=fread(fid, 1, 'int64', 'ieee-be')';

%         number of gridpoints along dimension 2: ndim2
ks.ndim2=fread(fid, 1, 'int64', 'ieee-be')';

%         number of gridpoints along dimension 3: ndim3
ks.ndim3=fread(fid, 1, 'int64', 'ieee-be')';

%     attribute mode: 8 byte integer (EAS3_NO_ATTR=1 or EAS3_ALL_ATTR=2)
EAS3_NO_ATTR=1; EAS3_ALL_ATTR=2;
ks.amode=fread(fid, 1, 'int64', 'ieee-be')';

%     geometry modes: each 8 byte integer (EAS3_NO_G=1, EAS3_X0DX_G=2, EAS3_UDEF_G=3, EAS3_ALL_G=4 or EAS3_FULL_G=5)
EAS3_NO_G=1; EAS3_X0DX_G=2; EAS3_UDEF_G=3; EAS3_ALL_G=4; EAS3_FULL_G=5;
%         time
ks.gmode(1)=fread(fid, 1, 'int64', 'ieee-be')';

%         parameter
ks.gmode(2)=fread(fid, 1, 'int64', 'ieee-be')';

%         spatial dimension 1
ks.gmode(3)=fread(fid, 1, 'int64', 'ieee-be')';

%         spatial dimension 2
ks.gmode(4)=fread(fid, 1, 'int64', 'ieee-be')';

%         spatial dimension 3
ks.gmode(5)=fread(fid, 1, 'int64', 'ieee-be')';


% array sizes for geometry data: each 8 byte integer
%     time
ks.ng(1)=fread(fid, 1, 'int64', 'ieee-be')';

%     parameter
ks.ng(2)=fread(fid, 1, 'int64', 'ieee-be')';

%     spatial dimension 1
ks.ng(3)=fread(fid, 1, 'int64', 'ieee-be')';

%     spatial dimension 2
ks.ng(4)=fread(fid, 1, 'int64', 'ieee-be')';

%     spatial dimension 3
ks.ng(5)=fread(fid, 1, 'int64', 'ieee-be')';

%     specification of user defined data: 8 byte integer (EAS3_NO_UDEF=1, EAS3_ALL_UDEF=2 or EAS3_INT_UDEF=3)
EAS3_NO_UDEF=1; EAS3_ALL_UDEF=2; EAS3_INT_UDEF=3;
ks.umode=fread(fid, 1, 'int64', 'ieee-be')';

%     array sizes for used defined data: each 8 byte integer
%         character field
ks.nu(1)=fread(fid, 1, 'int64', 'ieee-be')';

%         integer field
ks.nu(2)=fread(fid, 1, 'int64', 'ieee-be')';

%         real field
ks.nu(3)=fread(fid, 1, 'int64', 'ieee-be')';

%     time step array: nzs x 8 byte
ks.zsf=fread(fid, ks.nzs, 'int64', 'ieee-be')';

%     if attribute mode = EAS3_ALL_ATTR
if (ks.amode==EAS3_ALL_ATTR)
    %         time step attributes: nzs x ks.ATTRLEN x 1 byte character
    ks.zsattr=fread(fid, ks.nzs*ks.ATTRLEN, 'char*1', 'ieee-be')';
    ks.zsattr=char(reshape(ks.zsattr,ks.ATTRLEN,ks.nzs)');

    %         parameter attributes: npar x ks.ATTRLEN x 1 byte character
    ks.parattr=fread(fid, ks.npar*ks.ATTRLEN, 'char*1', 'ieee-be')';
    ks.parattr=char(reshape(ks.parattr,ks.ATTRLEN,ks.npar)');

    %         spatial attributes: 3 x ks.ATTRLEN x 1 byte character
    ks.dim1attr=fread(fid, ks.ATTRLEN, 'char*1', 'ieee-be')';
    ks.dim1attr=char(reshape(ks.dim1attr,ks.ATTRLEN,1)');

    ks.dim2attr=fread(fid, ks.ATTRLEN, 'char*1', 'ieee-be')';
    ks.dim2attr=char(reshape(ks.dim2attr,ks.ATTRLEN,1)');
    
    ks.dim3attr=fread(fid, ks.ATTRLEN, 'char*1', 'ieee-be')';
    ks.dim3attr=char(reshape(ks.dim3attr,ks.ATTRLEN,1)');
end

%     if geometry mode > EAS3_NO_G
%         for time, parameters and dimensions 1 to 3, each
%         geometry-array-size * 8 byte real
if (ks.gmode(1)>EAS3_NO_G)
    ks.gf(1).dat=fread(fid, ks.ng(1), 'float64', 'ieee-be')';
end    
if (ks.gmode(2)>EAS3_NO_G)
    ks.gf(2).dat=fread(fid, ks.ng(2), 'float64', 'ieee-be')';
end 
if (ks.gmode(3)>EAS3_NO_G)
    ks.gf(3).dat=fread(fid, ks.ng(3), 'float64', 'ieee-be')';
end
if (ks.gmode(4)>EAS3_NO_G)
    ks.gf(4).dat=fread(fid, ks.ng(4), 'float64', 'ieee-be')';
end
if (ks.gmode(5)>EAS3_NO_G)
    ks.gf(5).dat=fread(fid, ks.ng(5), 'float64', 'ieee-be')';
end

%     if user-defined data is chosen
%         character field: character-array-size x ks.UDEFLEN x 1 byte character
%         integer field: integer-array-size x 8 byte integer
%         real field: real-array-size x 8 byte real 
if (ks.umode==EAS3_ALL_UDEF)
    for i=1:ks.nu(1)
        ks.cf(i,:)=fread(fid, ks.UDEFLEN, 'uint8=>char', 'ieee-be')';
    end
    for i=1:ks.nu(2)
        ks.if(i)=fread(fid, 1, 'int64', 'ieee-be')';
    end
    for i=1:ks.nu(3)
        ks.rf(i)=fread(fid, 1, 'float64', 'ieee-be')';
    end    
end

data = zeros(ks.nzs, ks.npar, ks.ndim1, ks.ndim2, ks.ndim3);
tmp  = zeros(ks.ndim1*ks.ndim2*ks.ndim3,1);

% read data from file
for zs=1:ks.nzs
   for par=1:ks.npar
     tmp = fread(fid, ks.ndim1*ks.ndim2*ks.ndim3, 'float64', 'ieee-be');
     data(zs,par,:,:,:) = reshape(tmp, ks.ndim1, ks.ndim2, ks.ndim3);
   end
end    

fclose(fid);
