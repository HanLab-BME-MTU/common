classdef TestMultipleProgressText < TestCase
% Test the function multipleProgressText
    properties
    end
    methods
        function self = TestMultipleProgressText(name)
            self@TestCase(name);
        end
        function setUp(self)
            clear multipleProgressText;
        end
        function tearDown(self)
            clear multipleProgressText;
        end
        function testNormalCountdown(self)
            % Test a normal progress update scenario
            nSteps = [5 3 2];
            text = 'Testing';
            o = multipleProgressText(text,nSteps(1));
            assertEqual(o.level,1);
            assertEqual(o.iStep,0);
            assertEqual(o.nStep,nSteps(1));
            for i=1:nSteps(1)
                o = multipleProgressText(text,nSteps(2));
                assertEqual(o.level,2);
                assertEqual(o.iStep(2),0);
                assertEqual(o.nStep(2),nSteps(2));
                for j=1:nSteps(2)
                    o = multipleProgressText(text,nSteps(3));
                    assertEqual(o.level,3);
                    assertEqual(o.iStep(3),0);
                    assertEqual(o.nStep(3),nSteps(3));
                    for k=1:nSteps(3)
                        pause(0.1);
                        o = multipleProgressText;
                        assertEqual(o.iStep(3),k);
                        assertEqual(o.nStep(3),nSteps(3));
                    end
                    o = multipleProgressText;
                    assertEqual(o.iStep(2),j);
                    assertEqual(o.nStep(2),nSteps(2));
                end
                o = multipleProgressText;
                assertEqual(o.iStep(1),i);
                assertEqual(o.nStep(1),nSteps(1));
            end
        end
        function testOverrun(self)
            % Test situation where the number of steps specified is
            % exceeded
            nSteps = 5;
            overrun = 3;
            text = 'Testing';
            o = multipleProgressText(text,nSteps(1));
            assertEqual(o.level,1);
            assertEqual(o.iStep,0);
            assertEqual(o.nStep,nSteps(1));
            for i=1:(nSteps+overrun)
                pause(0.1);
                o = multipleProgressText;
            end
        end
    end
end