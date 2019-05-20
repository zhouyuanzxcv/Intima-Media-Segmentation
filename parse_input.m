function varargout = parse_input(varargin)
if ~isempty(varargin)
    if ischar(varargin{1})
        disp(['Process ',varargin{1}]);
        I = imread(varargin{1});
    elseif isa(varargin{1},'uint8')
        I = varargin{1};
    end
    varargout{1} = I;
end
if length(varargin) > 1
    option = varargin{2};
end
if isstruct(option)
    arg_set = fieldnames(option);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=option.',arg_set{i},';']);
        disp(['Use ',arg_set{i},' = ',num2str(option.(arg_set{i})(:)')]);
        eval(['varargout{',num2str(i+1),'}=',arg_set{i},';']);
    end
end
