classdef SOFEClass < handle
  properties
  end
  methods % constructor
    function obj = SOFEClass()
    end
  end
  methods
    function outputBig(obj, str)
      fprintf('#++++++++++++++++++++++++++\n');
      fprintf(['#  ', str, '\n']);
      fprintf('#++++++++++++++++++++++++++\n');
    end
  end
end