function igesout(obj,filename,varargin)
% IGESOUT IGES Converter for points, lines, and Nurbs curves and surfaces
%
% IGESOUT(OBJ,FILENAME) writes an IGES (*.igs) file with the given
% FILENAME from an object OBJ. OBJ can either be an Nx3 array of points,
% where each column corresponds to X, Y and Z respectively, or a
% structure describing a NURBS curve or surface.
%
% In addition OBJ may also be a cell array containing any combination of
% the above entities. IGESOUT will determine which type of IGES entity
% to encode based on the size and type of input. Therefore it is not
% necessary to specify what types of entities are being used.
%
% IGESOUT(OBJ,FILENAME,'PARAMETER','VALUE') allows the user to set a
% number of parameter/value pairs listed below
%
% PARAMETER DESCRIPTION
% --------- ------------------------------------------
% paramdelim Parameter Deliminator Character
% recorddelim Record Delimiter Character
% sendid Product ID from Sender
% receiveid Product ID for Receiver
% filename File Name
% systemid System ID
% preprocessor Pre-processor Version
% intbits Number of Bits for Integers
% precision1 Single Precision Magnitude
% significant1 Single Precision Significance
% precision2 Double Precision Magnitude
% significant2 Double Precision Significance
% scale Model Space Scale
% unitflag Unit Flag (1 = inches, 3 = look to units)
% units Model Units 'MM' or 'IN' (Default)
% lineweights Maximum number of line weights
% linewidth Maximum line width
% timestamp Time stamp of creation
% resolution Minimum User-inted Resolution
% max Approximate Maximum Coordinate
% author Author
% company Author's Organization
% version IGES Version Number
% standard Drafting Standard Code (- ANSI)
% date Model Creation/Change Date
%
% The defaults values for each of these parameters is already set within
% IGESOUT. It is not necessary to change them, as the default values
% seem to work fine. However, it may be prudent for some users to make
% modifications to these parameters. For instance, the default value for
% units is 'IN' (inches). To change to millimeters set units to 'MM'.
%
% It should also be noted that igesout also attempts to write the IGES
% file in the most compact form possible. This may make the file less
% readable, but it should minimize overall file size.
%
% In order to write Nurbs curves and surfaces to IGES, it will be
% neccessary to download and install the NURBS toolbox. It can be found
% at:
% http://www.mathworks.com/matlabcentral/fileexchange/
% <a href="http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=312&objectType=file">NURBS TOOLBOX</a>
%
%
% There are no guarantees that this program will work correctly. It was
% only tested with ProE and MasterCam, and seemed to work fine for each
% of the entities described above. If you have any problems, please send
% feedback.
%
% Example
% % Works if you have the NURBS toolbox (Saved in current directory)
% nrb = nrbtestsrf;
% crv = nrbtestcrv;
% pline = [0,0,1.1112;0,1.3333,2.1082;0,2.6667,1.6397;0,4,1.3346;];
% lin = pline(1,:)*3;
% pt = pline(1,:);
% igesout({nrb crv pline lin pt},'TestIGES_nurbs')
%
% % Works if you DON'T have the NURBS toolbox (Saved in current directory)
% pline = [0,0,1.1112;0,1.3333,2.1082;0,2.6667,1.6397;0,4,1.3346;];
% lin = pline(1,:)*3;
% pt = pline(1,:);
% igesout({pline lin pt},'TestIGES')

% Created by: Daniel Claxton
% dclaxton@ufl.edu
% 15-Mar-2007

% v1.01
% 08-Mar-2009
% Major speed improvements for large files

% Inspired by: Michael Fassbind's
% NrbsSrf2IGES.m
% 2006
%
% Ty - Edit 1/17/2021
% Updated to assumbe the nurbs objects have u,v 2 element vectors
% for the range of these parameters.

% Format input objects and determine Entity types
[obj,Entity,Num] = getEntities(obj);

% Decompose file name
[path, file]=fileparts(filename);
% if isempty(path)
% path=[cd,'\3D_parts'];
% end
% filename=fullfile(path, [file '.igs']);

% SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
% Start (S-SECTION)
% SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
S = ['Matlab to IGES converter. Written by Daniel Claxton',...
blanks(21) 'S0000001'];
% SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

% GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
% Setup Header (G-SECTION)
% GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
G.paramdelim = ','; % Parameter Deliminator Character
G.recorddelim = ';'; % Record Delimiter Character
G.sendid = 'Matlab'; % Product ID from Sender
G.receiveid = 'Receiver ID'; % Product ID for Receiver
G.filename = [file '.igs']; % File Name
G.systemid = computer; % System ID
G.preprocessor = 'Matlab -> IGES'; % Pre-processor Version
G.intbits = 16; % Number of Bits for Integers
G.precision1 = 06; % Single Precision Magnitude
G.significant1 = 15; % Single Precision Significance
G.precision2 = 13; % Double Precision Magnitude
G.significant2 = 15; % Double Precision Significance
G.scale = 1.0; % Model Space Scale
G.unitflag = 3; % Unit Flag (1 = inches, 3 = look to units)
G.units = 'IN'; % Model Units 'MM' or 'IN' (Default)
G.lineweights = 8; % Maximum number of line weights
G.linewidth = 0.016; % Maximum line width
G.timestamp = datestr(now); % Time stamp of creation
G.resolution = 1E-4; % Minimum User-inted Resolution
G.max = getMax(obj); % Approximate Maximum Coordinate
G.author = 'Author'; % Author
G.company = 'Authors Organization'; % Author's Organization
G.version = 11; % - USPRO/IPO-100 (IGES 5.2) [USPRO93]'; % IGES Version Number ** prob not right **
G.standard = 3; % - ANSI'; % Drafting Standard Code
G.date = datestr(now); % Model Creation/Change Date

% Check inputs
if nargin < 3
GG = G;
else
% Set above parameters from input parameter value pairs
GG = setParams(G,varargin{:});
end

G = G2str(GG);
% GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

%PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
% P-SECTION
%PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
P = zeros(14,Num);
Data = cell(Num,1);
for i = 1:Num;

if Entity(i) == 126 || Entity(i) == 128
X = squeeze(obj{i}.nurbs.coefs(1,:,:)); % Extract X Data
Y = squeeze(obj{i}.nurbs.coefs(2,:,:)); % Extract Y Data
Z = squeeze(obj{i}.nurbs.coefs(3,:,:)); % Extract Z Data
W = squeeze(obj{i}.nurbs.coefs(4,:,:)); % Extract Weights from coefficients

Coefs = [X(:) Y(:) Z(:)]';
Weights = W(:);
Coefs = Coefs(:);
end

switch Entity(i)
case 128
P(1,i) = Entity(i); % Entity Type. 128 = Rational spline surface
P(2,i) = obj{i}.k1; % Number of points in u direction (doesn't work without -1)
P(3,i) = obj{i}.k2; % Number of points in v direction
P(4,i) = obj{i}.m1; % Degree_u (degree = order -1)
P(5,i) = obj{i}.m2; % Degree_v
P(6,i) = obj{i}.prop1; % PROP1 Closed_u (0 = Not Closed)
P(7,i) = obj{i}.prop2; % PROP2 Closed_v (0 = Not Closed)
P(8,i) = obj{i}.prop3; % PROP3 (1 = Polynomial i.e. all weights equal, 0 = rational)***
P(9,i) = obj{i}.prop4; % PROP4 1st Direction periodicity (0 = Non-periodic)
P(10,i)= obj{i}.prop5; % PROP5 2nd Direction periodicity (0 = Non-periodic)
P(11,i)= obj{i}.u(1);
P(12,i)= obj{i}.u(2);
P(13,i)= obj{i}.v(1);
P(14,i)= obj{i}.v(2);
Uknots = obj{i}.s';
Vknots = obj{i}.t';
Knots = [Uknots; Vknots];
Data{i} = [Knots; Weights; obj{i}.q];
case 126
P(1,i) = Entity(i); % Entity Type. 126 = Rational spline curve
P(2,i) = obj{i}.number - 1; % Number of points - 1
P(3,i) = obj{i}.order - 1; % Degree of curve
P(4,i) = 0; % PROP1 Planar curve (0 = Not Planar)***
P(5,i) = 0; % PROP2 Closed curve (0 = Not Closed)
P(6,i) = 0; % PROP3 (1 = all weights equal, 0 = rational)***
P(7,i) = 0; % PROP4 Periodic (0 = Not periodic)
Knots = obj{i}.knots(:);
Data{i} = [Knots; Weights; Coefs];
case 106
P(1,i) = Entity(i); % Entity Type. 106 = Polyline
P(2,i) = 2; % Planar Polyline (2 = Not Planar)
P(3,i) = size(obj{i},1); % Number of points
% P(4,i) = 0; %
XYZ = obj{i}';
Data{i} = XYZ(:);
case 110
P(1,i) = Entity(i);
XYZ = obj{i}';
Data{i} = XYZ(:);
case 116
P(1,i) = Entity(i);
XYZ = obj{i}';
Data{i} = XYZ(:);
otherwise
end
% *** Note, we choose 0 because it is general case, it just requires us
% to write more information to the IGES file. But we can live with
% that to simplify things.

end
[Pstr,L] = P2str(P,Data,GG,Num);
%PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP

% DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
% D-SECTION
% DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

D = zeros(15,Num);
for i=1:Num
D(1,i) = Entity(i); % Entity Type.
D(2,i) = 1+sum(L(1:i)); % Data Start line
D(3,i) = 0; % Structure
D(4,i) = 1; % Line Font Pattern (1 = Solid)
D(5,i) = 0; % Level
D(6,i) = 0; % View
D(7,i) = 0; % Transformation Matrix
D(8,i) = 0; % Label Display
D(9,i) = 0; % Blank Status (0 = Visible)
D(10,i)= 0; % Subord. Entity Switch (0 = Independant)
D(11,i)= 0; % Entity Use Flag (0 = Geometry)
D(12,i)= 1; % Hierarchy ( 1 = Global defer)
D(13,i)= 0; % Line Weight Number
D(14,i)= L(i+1); % Data end line (Will be set later)
D(15,i)= 0; % Form Number (9 = General Quadratic Surface), 0 = none of above (1-9) options | *Modified Later
end
Dstr = D2str(D);
% DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

% 3DSKETCH 402
findpoint=find(Entity==116);
if ~isempty(findpoint)
if length(findpoint)==1
xa=1;
xb=1;
lx=1;
num=1;
else
[xa,xb,~,~,lx,num] = find_headtail(findpoint,1);
end
Nent=length(Entity);
sumL=sum(L);
nsp = 7;
NentD=D(2,end);
sumDL=size(Dstr,1);
for inum=1:num
Nent=Nent+1;
xc=xa(inum):xb(inum);

% xd=D(2,xc);

str1=sprintf('402,%1.0f,', lx(inum));
str2=sprintf('%1.0f,',2*xc-1);
str3=';';

Pstrtmp = [str1 str2 str3];

% Wrap String to be
Pstrtmp = wrap(Pstrtmp,GG.paramdelim,71-nsp);

[m] = size(Pstrtmp);

offset = sumL+(1:m);
Pset=offset*0+2*Nent-1;
count = sprintf('% 8.0fP%7.0f,',[Pset(:),offset(:)]');
count = wrap(count,',',17);

Pstrtmp = [Pstrtmp count(:,1:end-1)]; %#ok<AGROW>

Pstr=[Pstr; Pstrtmp];
L=[L; m];
sumL=sumL+m;

NentD=NentD+1;
sumDL=sumDL+2;
str1=sprintf(' 402% 8.0f 0 0 0 00000000D%7.0f', NentD,sumDL-1);
str2=sprintf(' 402 0 0% 8.0f 15 3DSKETCH 1D%7.0f', m,sumDL);

Dstr=[Dstr; str1;str2];
end
end

% TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
% T - SECTION
% TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
T = [size(S,1),size(G,1),size(Dstr,1),size(Pstr,1)];
T = sprintf('S%7.0fG%7.0fD%7.0fP%7.0f%sT%7.0f',T,blanks(40),1);
% TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

% Check for errors
chk = [size(S); size(G); size(Dstr); size(Pstr); size(T)];
if sum(diff(chk(:,2))) ~= 0
error(['Something went wrong writing the IGES file: ' file '.igs'])
end

% Concatenate all Sections
str = [S; G; Dstr; Pstr; T];

% Write to file
fid = fopen(filename,'w');
for i=1:size(str,1)-1
fprintf(fid,'%s\n',str(i,:));
end
fprintf(fid,'%s',str(i+1,:));
if fclose(fid) == 0
fprintf('%s\n IGS file written successfully!\n',filename);
end

function [xa,xb,diffa,diffb,lx,num] = find_headtail(x,n)
if isempty(x)
xa=[];
xb=[];
diffa=[];
diffb=[];
lx=0;
num=0;
return;
end

if nargin==1
n=1;
end

diffx=find(diff(x)>n);

[~, xd]=size(x);
if xd==1
diffa=[1;diffx+1];
diffb=[diffx;length(x)];
else
diffa=[1,diffx+1];
diffb=[diffx,length(x)];
end

xa=x(diffa);
xb=x(diffb);

lx=diffb-diffa+1;

num=length(xa);

function [obj,ent,n] = getEntities(obj)
if ~iscell(obj)
obj = {obj};
end
n = length(obj);
ent = zeros(1,n);
for i=1:n
ent(i)=obj{i}.type;
end

function G = setParams(G,p)
n = length(p);
n = n/2;
p = reshape(p,n,2);
for i=1:n
% G = setfield(lower(G,p{i,1}),p{i,2});
G = setfield(G,lower(p{i,1}),p{i,2});
end

function s = wrap(str,delim,n)
s = wrapmat(str,delim,n); return

function str = wrapmat(str,delim,nchar)
% Form the desired regular expression from nchars.
exp = sprintf('(.{1,%d})(?:%c+|$)', nchar-1, delim);

tokens = regexp(str,exp,'tokens');

% Each element if the cell array tokens is single-element cell array
% containing a string. Convert this to a cell array of strings.
get_contents = @(f) [f{1} delim];
c = cellfun(get_contents, tokens, 'UniformOutput', false);
c{end}(end) = [];

% Convert cell array into array of strings
str = cell2str(c);

% Insert extra spaces to assure strings are all [nchar] in length
nsrt = nchar-size(str,2);
if nsrt
% str(:,end+nsrt) = ' ';
str(:,end+1:end+nsrt)=' ';
end

function str = cell2str(c)
% Convert a cell array of strings into an array of strings.
% CELL2STR pads each string in order to force all strings
% have the same length.
%

% Determine the length of each string in cell array c
nblanks = cellfun(@length, c);
maxn = max(nblanks);
nblanks = maxn-nblanks;

% Create a cell array of blanks. Each column of the cell array contains
% the number of blanks necessary to pad each row of the converted string
padding = cellfun(@blanks,num2cell(nblanks), 'UniformOutput', false);

% Concatinate cell array and padding
str = {c{:}; padding{:}};

% This operation converts new the cell array into a string
str = [str{:}];

% Reshape the string into an array of strings
ncols = maxn;
nrows = length(str)/ncols;
str = reshape(str,ncols,nrows)';

function str = D2str(D)

k = 1;
[m,n]=size(D);
str = char(zeros(2*n,80-16));
for Nent=1:n
Entity = D(1,Nent);
switch Entity
case 128
case 126
case 110
case 116
case 106
D(15,Nent) = 12; % Change from Copious Data to Linear Curve
otherwise
end
tmp1 = sprintf('%8.0f', D(1:7,Nent));
tmp2 = sprintf('%8.0f',[D(1,Nent); D(12:15,Nent)]);
str(k,:) = [tmp1 blanks(23-8-7)];
k = k + 1;
str(k,:) = [tmp2 blanks(39-8-7)];
k = k + 1;
end

m = size(str,1);
on = true;
count = char(zeros(m,16));
for i=1:m
if on
count(i,:) = sprintf('% 8.0fD%7.0f',i,i);
else
count(i,:) = sprintf('% 8.0fD%7.0f',0,i);
end
on = ~on;
end

str = [str count];

function [str,L] = P2str(P,Data,G,Nentities)
pd = G.paramdelim;
rd = G.recorddelim;
prec = G.precision1;
str = [];
L = zeros(Nentities+1,1);
for Nent=1:Nentities
switch P(1,Nent)
case 128
plen = 1:10;
case 126
plen = 1:7;
case 106
plen = 1:3;
case 110
plen = 1;
case 116
plen = 1;
otherwise
end
str1 = sprintf(['%1.0f' pd],P(plen,Nent));
%precision = sprintf('%d',prec);
str2 = sprintf(['%1.6f' pd],Data{Nent});
str3 = sprintf(['%1.6f' rd],P(11:14,Nent));

% str3 = ['0., 2., 0., 2.' rd];

tmp = [str1 str2 str3];
nsp = 7;

% Wrap String to be
tmp = wrap(tmp,pd,71-nsp);

[m,n] = size(tmp);
L(Nent+1,1) = m;

for i=1:n-64
tmp(:,end+1) = ' '; %#ok<AGROW>
end

offset = sum(L(1:Nent));
for j=1:m
count(j,:) = sprintf('% 8.0fP%7.0f',2*Nent-1,j+offset);
end

tmp = [tmp count]; %#ok<AGROW>

str = [str; tmp];

count = []; % Clear count

end

function str = G2str(G)
fnames = fieldnames(G);
n = length(fnames);
str = cell(20,1);
pd = G.paramdelim;
rd = G.recorddelim;

k = 1;
for i=1:n
tmp = getfield(G,fnames{i}); %#ok<GFLD>
if ischar(tmp)
tmp = sprintf('%1.0fH%s%s',length(tmp),tmp,pd);
else
tmp = sprintf('%1.0f%s',tmp,pd);
end
if length(str{k}) + length(tmp) > 79-15
k = k + 1;
end

str{k,1} = [str{k} tmp];
if i == n
    str{k}(end) = rd;
end
end

last = 0;
while true
if isempty(str{last+1})
break
end
last = last+1;
end

str = str(1:last,1);
n = length(str{end,1});
stop = 80-8;
str{end}(n+1:stop) = blanks(stop-n);
str = char(str);

count = char(zeros(last,8));
for i=1:last
count(i,:) = sprintf('G%7.0f',i);
end

str = [str count];

function m = getMax(obj)
m = zeros(1,length(obj));
for i=1:length(obj)
if isstruct(obj{i})
m(i) = max(obj{i}.p(:));
else
m(i) = max(obj{i}(:));
end
end
m = max(m);