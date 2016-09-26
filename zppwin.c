// Part of Zinc FE package. Author: John Blackburn

#include <stdlib.h>
#include <windows.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <windowsx.h>

#include "zppwin.h" 

HINSTANCE g_hInstance;
char graphfile[1000];
int graphopen=0;
int iIndex,graphindex,ngraph;

char StemNameFull[MAX_PATH];
char command[MAX_PATH+4]="";

BOOL CALLBACK DlgProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam);
BOOL CALLBACK NewProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam);
int spliteq(char* str, char* left, char* right);
BOOL LoadTextFileToEdit(HWND hEdit, LPCTSTR pszFileName);
int getline(char s[], int lim, FILE *fp);
void Thread (PVOID pvoid);

typedef struct{
  HWND hwnd;
  // other stuff...
} PARAMS, *PPARAMS;

// ######################################################################

BOOL CALLBACK NewProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam)
{
  HDC hdc;
  HWND hCtrl;
  RECT rect;
  HENHMETAFILE hemf;
  
  if (Message==WM_INITDIALOG){
    graphopen=1;
    graphindex=iIndex;
    return TRUE;
  }

  else if (Message==WM_CLOSE){
    graphopen=0;
    EndDialog(hwnd, 0);
    return TRUE;
  }

  else if (Message==WM_COMMAND){
    if (LOWORD(wParam)==IDC_NEXT){
      //      printf("Hit Prev%d %d\n",graphindex,ngraph);
      if (graphindex<ngraph){
	graphindex++;
	sprintf(graphfile,"%s%02d.emf",StemNameFull,graphindex);
      }
    }
    else if (LOWORD(wParam)==IDC_PREV){
      if (graphindex>1){
	graphindex--;
	sprintf(graphfile,"%s%02d.emf",StemNameFull,graphindex);
      }
    }
    hCtrl=GetDlgItem(hwnd,IDC_PIC);
    GetClientRect(hCtrl,&rect);

    InvalidateRect (hCtrl, NULL, TRUE) ;
    UpdateWindow (hCtrl) ;
    
    hdc = GetDC (hCtrl) ;
    hemf=GetEnhMetaFile(graphfile);

    if (hemf != NULL){
    
      FillRect(hdc,&rect,GetStockObject(WHITE_BRUSH));
      PlayEnhMetaFile(hdc,hemf,&rect);
    
      DeleteEnhMetaFile(hemf);
    }
    else
      TextOut(hdc,rect.right/2,rect.bottom/2,"No file",7);

    TextOut(hdc,10,rect.bottom-20,graphfile,strlen(graphfile));
    ReleaseDC(hCtrl,hdc);


    return TRUE;
  }

  else if (Message==WM_PAINT){
    hCtrl=GetDlgItem(hwnd,IDC_PIC);
    GetClientRect(hCtrl,&rect);

    InvalidateRect (hCtrl, NULL, TRUE) ;
    UpdateWindow (hCtrl) ;
    
    hdc = GetDC (hCtrl) ;
    hemf=GetEnhMetaFile(graphfile);
    
    if (hemf !=NULL){
      FillRect(hdc,&rect,GetStockObject(WHITE_BRUSH));
      PlayEnhMetaFile(hdc,hemf,&rect);
    
      DeleteEnhMetaFile(hemf);
    }
    else
      TextOut(hdc,rect.right/2,rect.bottom/2,"No file",7);

    TextOut(hdc,10,rect.bottom-20,graphfile,strlen(graphfile));
    ReleaseDC(hCtrl,hdc);
    return FALSE;
  }

  else
    return FALSE;
}

// ######################################################################

BOOL CALLBACK DlgProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam)
{

  OPENFILENAME ofn;
  char PathName[MAX_PATH]="";
  char Path[MAX_PATH]="";

  char FileName[MAX_PATH]="";
  char OutName[MAX_PATH]="";
  char GnuName[MAX_PATH]="";
  char ZinName[MAX_PATH]="";
  char MeshName[MAX_PATH]="";
  char ConstName[MAX_PATH]="";
  char StemName[MAX_PATH]="";

  static char FileNameFull[MAX_PATH];
  static char OutNameFull[MAX_PATH];
  static char GnuNameFull[MAX_PATH];
  static char MeshNameFull[MAX_PATH];
  static char ZinNameFull[MAX_PATH];
  static char ConstNameFull[MAX_PATH];
  
  static HWND hOpen,hCancel,hRun,hVfile,hVgnu,hVlst,hEdit,
    hVmesh,hVconst,hVzin,hGraph;

  char string[1000],string2[1000];

  static PARAMS params;

  int i,plen,ist,gotfiles,rdu,ret,ind;
  FILE* fp;

  // ----------------------------------------------------------------------

  switch(Message){
  case WM_INITDIALOG:
    // This is where we set up the dialog box, and initialise any default values
    
    hRun=GetDlgItem(hwnd,IDC_RUN);     // all the buttons
    hOpen=GetDlgItem(hwnd,IDC_OPEN);     
    hCancel=GetDlgItem(hwnd,IDC_CANCEL);     
    hVfile=GetDlgItem(hwnd,IDC_VFILE);
    hVgnu=GetDlgItem(hwnd,IDC_VGNU);
    hVlst=GetDlgItem(hwnd,IDC_VLST);
    hVmesh=GetDlgItem(hwnd,IDC_VMESH);
    hVconst=GetDlgItem(hwnd,IDC_VCONST);
    hVzin=GetDlgItem(hwnd,IDC_VZIN);
    hEdit = GetDlgItem(hwnd, IDC_EDIT);  // the edit box
    hGraph=GetDlgItem(hwnd,IDC_GRAPH);

    SendMessage(hEdit,WM_SETFONT,(WPARAM)GetStockObject(SYSTEM_FIXED_FONT),TRUE);

    Button_Enable(hRun,FALSE);          // deactivate buttons and listbox
    Button_Enable(hVfile,FALSE);
    Button_Enable(hVgnu,FALSE);
    Button_Enable(hVlst,FALSE);
    Button_Enable(hVmesh,FALSE);
    Button_Enable(hVconst,FALSE);
    Button_Enable(hVzin,FALSE);
    Button_Enable(hGraph,FALSE);

    //    printf("%d",MAX_PATH);
    //    printf("size of OPENFILENAME %d",sizeof(OPENFILENAME));

    break;
  case WM_COMMAND:
    switch(LOWORD(wParam)){

      // ----------------------------------------------------------------------
      // Open .zin file and find the other files. View buttons initially greyed
      // ----------------------------------------------------------------------

    case IDC_OPEN:

      ZeroMemory(&ofn, sizeof(ofn));

      ofn.lStructSize = sizeof(OPENFILENAME);
      ofn.hwndOwner = hwnd;
      ofn.lpstrFilter = "Zinc Input Files (*.zpp)\0*.zpp\0";
      ofn.lpstrFile = PathName;
      ofn.nMaxFile = MAX_PATH;
      ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
      ofn.lpstrDefExt = "txt";

      // ----------------------------------------------------------------------
      // find the \ nearest to the end of the filename
      // Path is the string up to and including this slash
      // ----------------------------------------------------------------------

      if (GetOpenFileName(&ofn)){

	Button_Enable(hVfile,FALSE);
	Button_Enable(hVgnu,FALSE);
	Button_Enable(hVlst,FALSE);

	plen=strlen(PathName);

	ist=-1;
	for (i=plen;i>0;i--){
	  if (PathName[i]=='\\'){
	    ist=i;
	    break;
	  }
	}

	if (ist==-1) printf("could not find \\ \n");

	for (i=0;i<=ist;i++)
	  Path[i]=PathName[i];

	Path[ist+1]='\0';

	// ----------------------------------------------------------------------
	// FileName is the short name within the directory
	// MeshName is short name with .zin replaced with .mtf (ETC)
	// ----------------------------------------------------------------------

	strcpy(FileName,&PathName[ist+1]);

	strcpy(OutName,FileName);
	strcpy(GnuName,FileName);
	strcpy(MeshName,FileName);
	strcpy(ConstName,FileName);
	strcpy(ZinName,FileName);
	strcpy(StemName,FileName);

	for (i=strlen(FileName);i>0;i--){
	  if (strcmp(&FileName[i],".zpp")==0){
	    strcpy(&OutName[i],".zou");
	    strcpy(&GnuName[i],".gnu");
	    strcpy(&MeshName[i],".mtf");
	    strcpy(&ConstName[i],".con");
	    strcpy(&ZinName[i],".zin");
	    strcpy(&StemName[i],"");
	    break;
	  }
	}

	// ----------------------------------------------------------------------
	// FileNameFull (etc) includes the path
	// ----------------------------------------------------------------------

	strcpy(FileNameFull,Path);
	strcat(FileNameFull,FileName);
	//	printf("%s\n",FileNameFull);

	strcpy(OutNameFull,Path);
	strcat(OutNameFull,OutName);
	//printf("%s\n",OutNameFull);

	strcpy(GnuNameFull,Path);
	strcat(GnuNameFull,GnuName);
	//	printf("%s\n",GnuNameFull);

	strcpy(MeshNameFull,Path);
	strcat(MeshNameFull,MeshName);
	//printf("%s\n",MeshNameFull);

	strcpy(ConstNameFull,Path);
	strcat(ConstNameFull,ConstName);
	//printf("%s\n",ConstNameFull);

	strcpy(ZinNameFull,Path);
	strcat(ZinNameFull,ZinName);
	//printf("%s\n",ZinNameFull);

	strcpy(StemNameFull,Path);
	strcat(StemNameFull,StemName);
	//printf("%s\n",StemNameFull);

	// ----------------------------------------------------------------------
	// set "path" and "Filename" parts of control panel
	// ----------------------------------------------------------------------

	SetDlgItemText(hwnd, IDC_PATH, Path);
	SetDlgItemText(hwnd, IDC_FILENAME, FileName);

	// ----------------------------------------------------------------------
	// Go through the files and see if any are missing
	// ----------------------------------------------------------------------

	gotfiles=1;

	fp=fopen(FileNameFull,"r");
	if (fp==NULL){
	  SetDlgItemText(hwnd, IDC_FILENAME, "<File not found>");
	  gotfiles=0;
	}
	else {
	  fclose(fp);
	  SetDlgItemText(hwnd, IDC_FILENAME, FileName);
	  Button_Enable(hVfile,TRUE);
	}

	fp=fopen(OutNameFull,"r");
	if (fp==NULL){
	  SetDlgItemText(hwnd, IDC_OUTNAME, "<File not found>");
	  gotfiles=0;
	}
	else {
	  fclose(fp);
	  SetDlgItemText(hwnd, IDC_OUTNAME, OutName);
	}

	fp=fopen(MeshNameFull,"r");
	if (fp==NULL){
	  SetDlgItemText(hwnd, IDC_MESHNAME, "<File not found>");
	  gotfiles=0;
	}
	else {
	  fclose(fp);
	  SetDlgItemText(hwnd, IDC_MESHNAME, MeshName);
	  Button_Enable(hVmesh,TRUE);
	}

	fp=fopen(ConstNameFull,"r");
	if (fp==NULL){
	  SetDlgItemText(hwnd, IDC_CONSTNAME, "<No constants!>");
	}
	else {
	  fclose(fp);
	  SetDlgItemText(hwnd, IDC_CONSTNAME, ConstName);
	  Button_Enable(hVconst,TRUE);
	}

	fp=fopen(ZinNameFull,"r");
	if (fp==NULL){
	  SetDlgItemText(hwnd, IDC_ZINNAME, "<File not found>");
	  gotfiles=0;
	}
	else {
	  fclose(fp);
	  SetDlgItemText(hwnd, IDC_ZINNAME, ZinName);
	  Button_Enable(hVzin,TRUE);
	}

	// ----------------------------------------------------------------------
	// set gnuplot file text. 
	// ----------------------------------------------------------------------

	SetDlgItemText(hwnd, IDC_GNUNAME, GnuName);

	// ----------------------------------------------------------------------
	// Read file.zpp and set the LISTBOX
	// ----------------------------------------------------------------------

	SendDlgItemMessage(hwnd, IDC_VLST, LB_RESETCONTENT, 0, 0);

	fp=fopen(FileNameFull,"r");

	if (fp != NULL){

	  //	  printf("starting\n");
	  
	  ind=0;
	  while(getline(string,1000,fp) != -1){
	    //	    printf("%s\n",string);
	    if (strstr(string,"linescan")!=NULL || 
		strstr(string,"planescan")!=NULL || 
		strstr(string,"surfint")!=NULL ||
		strstr(string,"volint")!=NULL)
	      {
	      ind++;
	      sprintf(string2,"%s%02d.out << ",StemName,ind);
	      strcat(string2,string);
	      SendDlgItemMessage(hwnd, IDC_VLST, LB_ADDSTRING, 0, (LPARAM)string2);
	    }
	  }
	  ngraph=ind;

	  //	  Button_Enable(hVlst,TRUE);
	}

	fclose(fp);

	// ----------------------------------------------------------------------
	// light up RUN button if we got the files
	// ----------------------------------------------------------------------

	if (gotfiles==1) 
	  Button_Enable(hRun,TRUE);
	else
	  Button_Enable(hRun,FALSE);
      }
      break;
      
    case IDC_RUN:

      // ----------------------------------------------------------------------
      // Run zpp
      // ----------------------------------------------------------------------

      SetDlgItemText(hwnd, IDC_STATUS, "RUNNING");
      Button_Enable(hVlst,FALSE);
      Button_Enable(hVgnu,FALSE);
      Button_Enable(hRun,FALSE);
      Button_Enable(hOpen,FALSE);
      Button_Enable(hCancel,FALSE);

      strcpy(command,"zpp \"");
      strcat(command,StemNameFull);
      strcat(command,"\"");

      printf(command);
      //      system(command);

      params.hwnd=hwnd;
      _beginthread(Thread,0,&params);

      break;

    case IDC_VFILE:

      // ----------------------------------------------------------------------
      // Examine zpp file
      // ----------------------------------------------------------------------

      LoadTextFileToEdit(hEdit, FileNameFull);
      break;

    case IDC_VGNU:

      // ----------------------------------------------------------------------
      // Examine gnu file
      // ----------------------------------------------------------------------

      LoadTextFileToEdit(hEdit, GnuNameFull);
      break;

    case IDC_VCONST:

      // ----------------------------------------------------------------------
      // Examine const file
      // ----------------------------------------------------------------------

      LoadTextFileToEdit(hEdit, ConstNameFull);
      break;

    case IDC_VZIN:

      // ----------------------------------------------------------------------
      // Examine zin file
      // ----------------------------------------------------------------------

      LoadTextFileToEdit(hEdit, ZinNameFull);
      break;

    case IDC_VMESH:

      // ----------------------------------------------------------------------
      // Examine mesh file
      // ----------------------------------------------------------------------

      LoadTextFileToEdit(hEdit, MeshNameFull);
      break;

    case IDC_VLST:

      // ----------------------------------------------------------------------
      // List is clicked. See which item has been selected
      // ----------------------------------------------------------------------

      if (HIWORD(wParam)==1){
	iIndex=SendMessage(hVlst, LB_GETCURSEL, 0, 0)+1;
	sprintf(string,"%s%02d.out",StemNameFull,iIndex);
	sprintf(graphfile,"%s%02d.emf",StemNameFull,iIndex);
	//	printf("%s\n",graphfile);

	LoadTextFileToEdit(hEdit, string);
	if (GetFileAttributes(graphfile)!=0xFFFFFFFF) 
	  Button_Enable(hGraph,TRUE);
	else
	  Button_Enable(hGraph,FALSE);
      }
      break;
    case IDC_GRAPH:
      if (graphopen==0)
	DialogBox(g_hInstance, MAKEINTRESOURCE(IDD_NEW), NULL, NewProc);
      break;
    case IDC_CANCEL:
      EndDialog(hwnd, 0);
      break;
    default:
      return FALSE;
    }
    break;
  case WM_CLOSE:
    EndDialog(hwnd, 0);
    break;

  case WM_CALC_DONE:
      SetDlgItemText(hwnd, IDC_STATUS, "READY");
      Button_Enable(hVlst,TRUE);
      Button_Enable(hVgnu,TRUE);
      Button_Enable(hRun,TRUE);
      Button_Enable(hOpen,TRUE);
      Button_Enable(hCancel,TRUE);
      break;
  default:
    return FALSE;
  }
  return TRUE;
}

// ######################################################################

int spliteq(char* str, char* left, char* right)
{

  // ----------------------------------------------------------------------
  // Split string str at = sign to get null-terminated left and right
  // remove leading blanks from both returned strings
  // ----------------------------------------------------------------------

  int i,j,ind;

  for (i=0;i<strlen(str);i++){
    if (str[i]=='='){
      ind=0;
      for (j=0;j<i;j++){
	if (str[j]!=' '){
	  left[ind]=str[j];
	  ind=ind+1;
	}
      }
      left[ind]='\0';
      
      ind=0;
      for (j=i+1;j<strlen(str);j++){
	if (str[j]!=' ' && str[j]!='\n'){
	  right[ind]=str[j];
	  ind=ind+1;
	}
      }
      right[ind]='\0';
      
      return 1;
    }
  }
  return 0;
}

// ######################################################################

BOOL LoadTextFileToEdit(HWND hEdit, LPCTSTR pszFileName)
{
  HANDLE hFile;
  BOOL bSuccess = FALSE;

  hFile = CreateFile(pszFileName, GENERIC_READ, FILE_SHARE_READ, NULL,
		     OPEN_EXISTING, 0, NULL);
  if(hFile != INVALID_HANDLE_VALUE)
    {
      DWORD dwFileSize;
      
      dwFileSize = GetFileSize(hFile, NULL);
      if(dwFileSize != 0xFFFFFFFF)
	{
	  LPSTR pszFileText;
	  
	  pszFileText = (LPSTR)GlobalAlloc(GPTR, dwFileSize + 1);
	  if(pszFileText != NULL)
	    {
	      DWORD dwRead;
	      
	      if(ReadFile(hFile, pszFileText, dwFileSize, &dwRead, NULL))
		{
		  pszFileText[dwFileSize] = 0; // Add null terminator
		  if(SetWindowText(hEdit, pszFileText))
		    bSuccess = TRUE; // It worked!
		}
	      GlobalFree(pszFileText);
	    }
	}
      CloseHandle(hFile);
    }
  return bSuccess;
}

// ######################################################################

int getline(char s[], int lim, FILE *fp)
{

  int and,i;
  char c;

  and=0;

  i=0;
  while (i<lim){
    c=fgetc(fp);
    if (c=='&') and=1;
    else if (c==EOF) break;
    else if (c=='\n'){
      if (and==0)
	break;
      else
	and=0;
    }
    else if (and==0){
      s[i]=tolower(c);
      i++;
    }
  }
  
  s[i]='\0';

  if (c==EOF)
    return -1;
  else
    return i;
  
}

void Thread (PVOID pvoid)
{

  PPARAMS pparams;

  pparams=(PPARAMS) pvoid;

  system(command);

  SendMessage(pparams->hwnd, WM_CALC_DONE, 0, 0);

  _endthread();
}

// ######################################################################

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
		   LPSTR lpCmdLine, int nCmdShow)
{
  g_hInstance=hInstance;
  return DialogBox(hInstance, MAKEINTRESOURCE(IDD_MAIN), NULL, DlgProc);
}

int main()
{
  return WinMain(GetModuleHandle(NULL),NULL,GetCommandLineA(),SW_SHOWNORMAL);
}
