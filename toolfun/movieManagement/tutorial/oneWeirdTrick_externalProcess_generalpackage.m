MD = MovieData.load('example.tif');
threshProc = ThresholdProcess(MD);
extProc = ExternalProcess(MD,'Say something',@(p) disp(p.getParameters().text));
extProc.setParameters(struct('text','Hello world!'));
MD.addProcess(threshProc);
MD.addProcess(extProc);
MD.addPackage(GenericPackage(MD));
MD.packages_{1}.GUI(MD);
