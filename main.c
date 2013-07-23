#include <stdlib.h>
#include <stdio.h>
#include <SDL/SDL.h>
#include <math.h>
#include "portaudio.h"

//# define TASK_NOTES {40,41,43,45,47,48} // C maj big
//# define TASK_NOTES {48,50,52,53,55,57,59,60} // C maj small
//# define TASK_NOTES {60,62,64,65,67,69,71,72} // C maj 1
//# define TASK_NOTES {72,74,76} // C maj 2
# define TASK_NOTES {40,41,43,45,47,48, 48,50,52,53,55,57,59,60, 62,64,65,67,69,71,72,74,76}


//for portaudio
#define SAMPLE_RATE (44100)
//       N    = 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384;
//       LogN = 2, 3,  4,  5,  6,   7,   8,   9,   10,   11,   12,   13,    14;
#define FRAMES_PER_BUFFER (8192)
#define POW_2 (13)
#define NUM_CHANNELS (1)
#define DELTA_NOTE 10

/* Select sample format. */
#define PA_SAMPLE_TYPE paFloat32
#define SAMPLE_SIZE (4)
#define SAMPLE_SILENCE (0.0f)
#define CLEAR(a) memset( (a), 0, FRAMES_PER_BUFFER * NUM_CHANNELS * SAMPLE_SIZE ) // обнуление массива
#define PRINTF_S_FORMAT "%.8f"

#define  NUMBER_IS_2_POW_K(x)   ((!((x)&((x)-1)))&&((x)>1))  // x is pow(2, k), k=1,2, ...
#define  FT_DIRECT        -1    // Direct transform.
#define  FT_INVERSE        1    // Inverse transform.

SDL_Surface *LINES;
SDL_Surface *LINEADD;
SDL_Surface *NOTE;
SDL_Surface *SCREEN;
SDL_Rect dest;
SDL_Rect destline;

float X0, Y0;

struct StructNote_
{
    int num;
    float freq;
    float freq_down;
    float freq_up;
} note[128];

void VerPar(float x1,float y1,float x2,float y2,float x3,float y3)
{
    // определяем координаты вершины параболы, выводим в глобал переменные
    float A,B,C;


    A=(y3-( (x3*(y2-y1)+x2*y1-x1*y2) /  (x2-x1)))   /(x3*(x3-x1-x2)+x1*x2);
    B= (y2-y1)/(x2-x1) - A*(x1+x2);
    C= ( (x2*y1-x1*y2) / (x2-x1) )+A*x1*x2;
    X0= -(B/(2*A));
    Y0= - ( (B*B-4*A*C)  /   (4*A) );
    //printf("\ny1=%g y2=%g y3=%g => Y0=%g \n",y1,y2,y3,Y0);
}

int  FFT(float *Rdat, float *Idat, int N, int LogN, int Ft_Flag)
{
    // parameters error check:
    if((Rdat == NULL) || (Idat == NULL))                  return 0;
    if((N > 16384) || (N < 1))                            return 0;
    if(!NUMBER_IS_2_POW_K(N))                             return 0;
    if((LogN < 2) || (LogN > 14))                         return 0;
    if((Ft_Flag != FT_DIRECT) && (Ft_Flag != FT_INVERSE)) return 0;

    register int  i, j, n, k, io, ie, in, nn;
    float         ru, iu, rtp, itp, rtq, itq, rw, iw, sr;

    static const float Rcoef[14] =
    {
        -1.0000000000000000F,  0.0000000000000000F,  0.7071067811865475F,
        0.9238795325112867F,  0.9807852804032304F,  0.9951847266721969F,
        0.9987954562051724F,  0.9996988186962042F,  0.9999247018391445F,
        0.9999811752826011F,  0.9999952938095761F,  0.9999988234517018F,
        0.9999997058628822F,  0.9999999264657178F
    };
    static const float Icoef[14] =
    {
        0.0000000000000000F, -1.0000000000000000F, -0.7071067811865474F,
        -0.3826834323650897F, -0.1950903220161282F, -0.0980171403295606F,
        -0.0490676743274180F, -0.0245412285229122F, -0.0122715382857199F,
        -0.0061358846491544F, -0.0030679567629659F, -0.0015339801862847F,
        -0.0007669903187427F, -0.0003834951875714F
    };

    nn = N >> 1;
    ie = N;
    for(n=1; n<=LogN; n++)
    {
        rw = Rcoef[LogN - n];
        iw = Icoef[LogN - n];
        if(Ft_Flag == FT_INVERSE) iw = -iw;
        in = ie >> 1;
        ru = 1.0F;
        iu = 0.0F;
        for(j=0; j<in; j++)
        {
            for(i=j; i<N; i+=ie)
            {
                io       = i + in;
                rtp      = Rdat[i]  + Rdat[io];
                itp      = Idat[i]  + Idat[io];
                rtq      = Rdat[i]  - Rdat[io];
                itq      = Idat[i]  - Idat[io];
                Rdat[io] = rtq * ru - itq * iu;
                Idat[io] = itq * ru + rtq * iu;
                Rdat[i]  = rtp;
                Idat[i]  = itp;
            }

            sr = ru;
            ru = ru * rw - iu * iw;
            iu = iu * rw + sr * iw;
        }

        ie >>= 1;
    }

    for(j=i=1; i<N; i++)
    {
        if(i < j)
        {
            io       = i - 1;
            in       = j - 1;
            rtp      = Rdat[in];
            itp      = Idat[in];
            Rdat[in] = Rdat[io];
            Idat[in] = Idat[io];
            Rdat[io] = rtp;
            Idat[io] = itp;
        }

        k = nn;

        while(k < j)
        {
            j   = j - k;
            k >>= 1;
        }

        j = j + k;
    }

    if(Ft_Flag == FT_DIRECT) return 1;

    rw = 1.0F / N;

    for(i=0; i<N; i++)
    {
        Rdat[i] *= rw;
        Idat[i] *= rw;
    }

    return 1;
}

void FreqOfNotes()
{
    int i;
    for (i=0; i<=127; i++)
    {
        note[i].num=i;
        note[i].freq = 440 * pow(2,(i-69)/12.0) ;
        note[i].freq_up  =note[i].freq* pow( pow (2,1/1200.), DELTA_NOTE);
        note[i].freq_down=note[i].freq* pow( pow (2,1/1200.), -DELTA_NOTE);
    }
    return ;
}

void DrawAddLines(int AddLines)
{
    destline.x = 340;

    switch (AddLines)
    {
    case -3:
        destline.y = 347;
        SDL_BlitSurface(LINEADD, NULL, SCREEN, &destline);
    case -2:
        destline.y = 317;
        SDL_BlitSurface(LINEADD, NULL, SCREEN, &destline);
    case -1:
        destline.y = 287;
        SDL_BlitSurface(LINEADD, NULL, SCREEN, &destline);
        break;
    case 3:
        destline.y = 20;
        SDL_BlitSurface(LINEADD, NULL, SCREEN, &destline);
    case 2:
        destline.y = 50;
        SDL_BlitSurface(LINEADD, NULL, SCREEN, &destline);
    case 1:
        destline.y = 80;
        SDL_BlitSurface(LINEADD, NULL, SCREEN, &destline);
    }
    return ;
}

void DrawNote(int CurNote)
{
    SDL_BlitSurface(LINES, NULL, SCREEN, NULL);

    dest.x = 360;
    switch (CurNote)
    {
        // 40,41,43,45,47,48  C maj big
    case 40:
        DrawAddLines(-3);
        dest.y = 350;
        break;
    case 41:
        DrawAddLines(-3);
        dest.y = 334;
        break;
    case 43:
        DrawAddLines(-2);
        dest.y = 320;
        break;
    case 45:
        DrawAddLines(-2);
        dest.y = 304;
        break;
    case 47:
        DrawAddLines(-1);
        dest.y = 290;
        break;

        //48,50,52,53,55,57,59,60  C maj small
    case 48:
        DrawAddLines(-1);
        dest.y = 274;
        break;
    case 50:
        dest.y = 258;
        break;
    case 52:
        dest.y = 241;
        break;
    case 53:
        dest.y = 222;
        break;
    case 55:
        dest.y = 204;
        break;
    case 57:
        dest.y = 186;
        break;
    case 59:
        dest.y = 168;
        break;
        // 60,62,64,65,67,69,71,72 C maj 1
    case 60:
        dest.y = 151;
        break;
    case 62:
        dest.y = 134;
        break;
    case 64:
        dest.y = 115;
        break;
    case 65:
        dest.y = 98;
        break;
    case 67:
        dest.y = 80;
        break;
    case 69:
        DrawAddLines(1);
        dest.y = 66;
        break;
    case 71:
        DrawAddLines(1);
        dest.y = 48;
        break;

        // 72,74,76 // C maj 2
    case 72:
        DrawAddLines(2);
        dest.y = 35;
        break;
    case 74:
        DrawAddLines(2);
        dest.y = 18;
        break;
    case 76:
        DrawAddLines(3);
        dest.y = 5;
        break;
    }

    SDL_BlitSurface(NOTE, NULL, SCREEN, &dest);
    SDL_Flip(SCREEN);
    return ;
}

int GetNote (float *Rdat, float *Idat, int N)
{
    float max_ampl, ampl_2, ampl_3;
    int i, max_i, i_2, i_3, note_real;
    int i_start=60*FRAMES_PER_BUFFER/SAMPLE_RATE;// нижняя частота
    int i_last=700*FRAMES_PER_BUFFER/SAMPLE_RATE;// верхняя частота
    float len_vector[FRAMES_PER_BUFFER];
    float Freq_1=0, Freq_2=0, Freq_3=0, FreqReal=0;
    // определяем частоту с максимальным вектором
            max_ampl=0;
            max_i=0;
            for(i=i_start; i<i_last; i++)
            {
                len_vector[i]=Rdat[i]*Rdat[i]+Idat[i]*Idat[i];
                if (len_vector[i]>max_ampl)
                {
                    max_ampl=len_vector[i];
                    max_i=i;
                }
            }
            // интерполяция, определяем вершину параболы для истиной частоты
            VerPar((float)max_i-1,
                   sqrt(  ( Rdat[max_i-1]*Rdat[max_i-1])+(Idat[max_i-1]*Idat[max_i-1]) ),
                   (float)max_i,
                   sqrt(  max_ampl ),
                   (float)max_i+1,
                   sqrt(  (Rdat[max_i+1]*Rdat[max_i+1])+(Idat[max_i+1]*Idat[max_i+1])));

            Freq_1=X0*SAMPLE_RATE/FRAMES_PER_BUFFER;
            Freq_2=Freq_1/2;
            Freq_3=Freq_1/3;

            // ищем индексы и мощности гармонических частот
            i_2=floor((X0/2)+0.5);
            ampl_2=sqrt(  ( Rdat[i_2]*Rdat[i_2])+(Idat[i_2]*Idat[i_2]) );
            //printf("Pow_2 %4.2f \t", ampl_2 );

            i_3=floor((X0/3)+0.5);
            ampl_3=sqrt(  ( Rdat[i_3]*Rdat[i_3])+(Idat[i_3]*Idat[i_3]) );
            //printf("Pow_3 %4.2f \n", ampl_3);

            // анализиурем по гармоникам
            FreqReal=Freq_1;
            if (ampl_2 > 60 && Freq_2>80) FreqReal=Freq_2;
            if (ampl_3 > 60 && Freq_3>80) FreqReal=Freq_3;

            // get number of note
            note_real=0;
            for (i=40; i<=83; i++)
            {
                if (FreqReal>note[i].freq_down && FreqReal<note[i].freq_up) note_real=i;
            }
    return (note_real);
}

int main ( int argc, char** argv )
{
    int CountTaskNote,CurNote,LastNote=0;
    int TaskNote[]=TASK_NOTES;
    int Solve=1;
    CountTaskNote= sizeof(TaskNote)/4;

    // initialize SDL video
    SDL_Init(SDL_INIT_VIDEO);
    // make sure SDL cleans up before exit
    atexit(SDL_Quit);
    // create a new window
    SCREEN = SDL_SetVideoMode(600,400,32, SDL_HWSURFACE|SDL_DOUBLEBUF);
    SDL_WM_SetCaption("C maj","ex1");

    // load images
    LINES=SDL_LoadBMP("NoteLines_600_400.bmp");
    NOTE=SDL_LoadBMP("Note.bmp");
    LINEADD=SDL_LoadBMP("LineAdd.bmp");

    //инициализация портаудио
    FreqOfNotes();// заполняем структуру нот
    PaStreamParameters inputParameters;
    PaStream *stream = NULL;
    float *sampleBlock;
    int j,Real_Note=0;
    int numBytes;
    float Im[FRAMES_PER_BUFFER];

    numBytes = FRAMES_PER_BUFFER * NUM_CHANNELS * SAMPLE_SIZE ;
    sampleBlock = (float *) malloc( numBytes );

    CLEAR( sampleBlock );

    Pa_Initialize();
    inputParameters.device = Pa_GetDefaultInputDevice();
    inputParameters.channelCount = NUM_CHANNELS;
    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultHighInputLatency ;
    inputParameters.hostApiSpecificStreamInfo = NULL;

    Pa_OpenStream(
        &stream,
        &inputParameters,
        0,
        SAMPLE_RATE,
        FRAMES_PER_BUFFER,
        paClipOff,
        NULL,
        NULL );

    Pa_StartStream( stream );
    fflush(stdout);

    // ##########################################
    // program main loop
    int done=0;
    while(done==0)
    {
        // message processing loop
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            // check for messages
            switch (event.type)
            {
            case SDL_QUIT:
                done = 1;
                break;
            case SDL_KEYDOWN:
            {
                if (event.key.keysym.sym == SDLK_ESCAPE)
                    done = 1;
                break;
            }

            }// end switch
        }// end of message processing

        // ##########################################
        // GENERATE AND DRAWING

        if (Solve==1)
        {
            do
            {
                CurNote=( rand()%(CountTaskNote));
            }
            while (CurNote==LastNote);
            LastNote=CurNote;
            DrawNote(TaskNote[CurNote]); //DRAW SDL
            Solve=0;
        }

        // ######################################
        // PortAudio analyzer

        float max_sample=0;
        Pa_ReadStream( stream, sampleBlock, FRAMES_PER_BUFFER );

        // определяем максимальную амплитуду
        for (j=0; j<FRAMES_PER_BUFFER; ++j)
        {
            if (sampleBlock[j] > max_sample) max_sample=sampleBlock[j];
        }

        // ***************************************
        // если аплтитуда мала, ничего не анализируем!
        if (max_sample>0.2)
        {
            // отправляем массив на БПФ, обнулив мнимую часть
            CLEAR( Im );
            FFT(sampleBlock, Im, FRAMES_PER_BUFFER, POW_2, -1);

            Real_Note=GetNote(sampleBlock, Im, FRAMES_PER_BUFFER);


        }
        if (Real_Note==TaskNote[CurNote]) Solve=1 ;

    }// end main loop

//close portaudio
    CLEAR( sampleBlock );
    free( sampleBlock );

    Pa_StopStream( stream );
    Pa_Terminate();

// free loaded bitmap
    SDL_FreeSurface(LINES);
    SDL_FreeSurface(LINEADD);
    SDL_FreeSurface(NOTE);

    return 0;
}
