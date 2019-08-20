//
//  main.cpp
//  ultraP
//


#include <iostream>
#include "FASTAReader.hpp"
#include "ultra.hpp"
#include "settings.hpp"

int main(int argc, const char * argv[]) {
    // insert code here...
    
    
   Settings settings = Settings(argc, argv);
  //  const char * argu[] ={"./ultra", "-m", "0.85", "-mi", "0.01", "-md", "0.01",  "/Users/daniel/Desktop/chr18.fa"};
   // Settings settings = Settings(8, argu);
    
   // settings.v_windowSize = 20000;
   // settings.v_maxInsertion = 10;
   // settings.v_maxDeletion = 5;
    
  //  printf("Reading ultra settings.\n");
    Ultra *ultra = new Ultra(&settings, 0);
  //  printf("Analyzing file.\n");
    ultra->AnalyzeFile();
  //  printf("Outputing remaining repeats.\n");
    ultra->OutputRepeats(true);
    
    if (ultra->settings->v_outFilePath != "")
        fclose(ultra->out);
  //  printf("Terminating.\n");
    //printf("%i\n", ultra->reader->doneReadingFile);
    
    
  /*  FASTAReader *reader = new FASTAReader(settings.v_filePath,
                                          settings.v_numberOfWindows,
                                          settings.v_windowSize,
                                          settings.v_overlapSize);
    
    
    int leng = settings.v_windowSize + settings.v_overlapSize + 2;
    
    UModel *model = new UModel(settings.v_maximumPeriod,
                               settings.v_maxInsertion,
                               settings.v_maxDeletion,
                               leng);
    
    model->SetMatchProbabilities(settings.v_matchProbability);
    
    model->tp_zeroToMatch = settings.v_zeroToMatch;
    model->tp_matchToZero = settings.v_matchToZero;
    
    model->tp_matchToInsertion  = settings.v_matchToInsertion;
    model->tp_matchToDeletion   = settings.v_matchToDeletion;
    
    model->tp_consecutiveInsertion  = settings.v_consecutiveInsertion;
    model->tp_consecutiveDeletion   = settings.v_consecutiveDeletion;
    
    model->CalculateScores();
    
    Ultra *ultra = new Ultra(model, reader);
    ultra->scoreThreshold = settings.v_scoreThreshold;
    ultra->outputReadID = true;
    ultra->outputRepeatSequence = true;
    ultra->AnalyzeFile();
    ultra->OutputRepeats();*/
    
    
    
  /*  for (int i = 0; i < ultra->repeats.size(); ++i) {
        RepeatRegion *r = ultra->repeats[i];
        printf("%s: %llu %llu: %f: \n", r->sequenceName.c_str(), r->sequenceStart, r->readID, r->regionScore);
        
        for (int j = 0; j < r->repeatPeriod; ++j) {
            printf("%c", CharForSymbol(r->consensus[j]));
        }
        
        printf("\n");
    }
    
    delete reader;
    
    
    
    
    SequenceWindow *window = reader->readyWindows.top();
    printf("%s\n", window->sequenceName.c_str());
    
  //  int cr = 0;
    while (reader->readyWindows.size() > 0) {
        SequenceWindow *window = reader->readyWindows.top();
        reader->readyWindows.pop();
        
        
        printf("\n\n-------------------------------------\n");
        printf("Reading window with...\n");
        printf("Sequence name: %s\n", window->sequenceName.c_str());
        printf("Sequence ID: %lli\n", window->seqID);
        printf("Read ID: %lli\n", window->readID);
        printf("Window ID: %lli\n", window->windowID);
        printf("Start: %lli end: %lli\n", window->start, window->end);
        
        printf("\n------------\nOverlap\n------------\n");
        
        for (int i = 0; i < window->overlap; ++i) {
            printf("%c", CharForSymbol(window->overlapSeq[i]));
        }
        printf("\n------------\nNew Seq\n------------\n");
        
        for (int i = 0; i < window->length; ++i) {
            printf("%c", CharForSymbol(window->newSeq[i]));
        }
        
        printf("\n------------\nWhole Seq\n------------\n");
        
        for (int i = 0; i < window->length + window->overlap; ++i) {
            printf("%c", CharForSymbol(window->seq[i]));
        }
    }
    */
    return 0;
}
