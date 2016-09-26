// Part of Zinc FE package. Author: John Blackburn

#include <stdlib.h>
#include <windows.h>
#include <string.h>
#include <stdio.h>
#include <windowsx.h>

#include "zincfe.h" 

BOOL CALLBACK DlgProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam);
int spliteq(char* str, char* left, char* right);
BOOL LoadTextFileToEdit(HWND hEdit, LPCTSTR pszFileName);
void Thread (PVOID pvoid);

char command[MAX_PATH+4]="";

int WM_calc_done;        // value set by RegisterWindowMessage

typedef struct{
  HWND hwnd;
  // add other data...
} PARAMS, *PPARAMS;

//----------------------------------------------------------------------

BOOL CALLBACK DlgProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam)
{

  static PARAMS params;

  OPENFILENAME ofn;
  char PathName[MAX_PATH]="";
  char Path[MAX_PATH]="";

  char FileName[MAX_PATH]="";
  char MeshName[MAX_PATH]="";
  char RstName[MAX_PATH]="";
  char OutName[MAX_PATH]="";
  char LstName[MAX_PATH]="";
  char ConstName[MAX_PATH]="";
  char NLName[MAX_PATH]="";

  static char FileNameFull[MAX_PATH];
  static char FileNameStem[MAX_PATH];
  static char MeshNameFull[MAX_PATH];
  static char RstNameFull[MAX_PATH];
  static char OutNameFull[MAX_PATH];
  static char LstNameFull[MAX_PATH];
  static char ConstNameFull[MAX_PATH];
  static char NLNameFull[MAX_PATH];
  static char extension[10]=".f90";

  static HWND hOpen,hCancel,hRun,hVfile,hVmesh,hVlst,hEdit,hVconst,hVNL;

  char string[100],left[100],right[100];

  int i,plen,ist,gotfiles,rdu,ret;
  FILE* fp;

  // ----------------------------------------------------------------------

  if (Message==WM_INITDIALOG){
    // This is where we set up the dialog box, and initialise any default values
    
    WM_calc_done=RegisterWindowMessage("zincfe_calc_done");   // this is not a constant!
    
    GetCurrentDirectory(MAX_PATH,PathName);

    printf("PWD=%s\n",PathName);

    //    printf("%d %d %d %d\n",WM_INITDIALOG,WM_USER,WM_COMMAND,WM_calc_done);
    
    SetDlgItemText(hwnd, IDC_EXT, extension);
    
    hOpen=GetDlgItem(hwnd,IDC_OPEN);
    hCancel=GetDlgItem(hwnd,IDC_CANCEL);
    hRun=GetDlgItem(hwnd,IDC_RUN);
    hVfile=GetDlgItem(hwnd,IDC_VFILE);
    hVmesh=GetDlgItem(hwnd,IDC_VMESH);
    hVlst=GetDlgItem(hwnd,IDC_VLST);
    hVconst=GetDlgItem(hwnd,IDC_VCONST);
    hVNL=GetDlgItem(hwnd,IDC_VNL);
    hEdit = GetDlgItem(hwnd, IDC_EDIT);
    
    SendMessage(hEdit,WM_SETFONT,(WPARAM)GetStockObject(SYSTEM_FIXED_FONT),TRUE);
    
    Button_Enable(hRun,FALSE);
    Button_Enable(hVfile,FALSE);
    Button_Enable(hVmesh,FALSE);
    Button_Enable(hVlst,FALSE);
    Button_Enable(hVconst,FALSE);
    Button_Enable(hVNL,FALSE);

      //    printf("%d",MAX_PATH);
      //printf("size of OPENFILENAME %d",sizeof(OPENFILENAME));
  }
  else if (Message==WM_COMMAND){
    
    // ----------------------------------------------------------------------
    // Open .zin file and find the other files. View buttons initially greyed
    // ----------------------------------------------------------------------
    
    if (LOWORD(wParam)==IDC_OPEN){
      
      Button_Enable(hVfile,FALSE);
      Button_Enable(hVmesh,FALSE);
      Button_Enable(hVlst,FALSE);
      Button_Enable(hVconst,FALSE);
      Button_Enable(hVNL,FALSE);
      
      ZeroMemory(&ofn, sizeof(ofn));
      
      ofn.lStructSize = sizeof(OPENFILENAME);
      ofn.hwndOwner = hwnd;
      ofn.lpstrFilter = "Zinc Input Files (*.zin)\0*.zin\0";
      ofn.lpstrFile = PathName;
      ofn.nMaxFile = MAX_PATH;
      ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
      ofn.lpstrDefExt = "txt";
	
      // ----------------------------------------------------------------------
      // find the \ nearest to the end of the filename
      // Path is the string up to and including this slash
      // ----------------------------------------------------------------------
      
      if (GetOpenFileName(&ofn)){
	
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
	  
	strcpy(MeshName,FileName);
	strcpy(RstName,FileName);
	strcpy(OutName,FileName);
	strcpy(LstName,FileName);
	strcpy(ConstName,FileName);
	strcpy(NLName,FileName);
	
	GetDlgItemText(hwnd, IDC_EXT, extension, 10);
	
	for (i=strlen(FileName);i>0;i--){
	  if (strcmp(&FileName[i],".zin")==0){
	    strcpy(&MeshName[i],".mtf");
	    strcpy(&RstName[i],".rst");
	    strcpy(&OutName[i],".zou");
	    strcpy(&LstName[i],".zls");
	    strcpy(&ConstName[i],".con");
	    strcpy(&NLName[i],extension);
	    break;
	  }
	}

	// ----------------------------------------------------------------------
	// FileNameFull (etc) includes the path
	// ----------------------------------------------------------------------
	  
	strcpy(FileNameFull,Path);
	strcat(FileNameFull,FileName);
	//	printf("%s\n",FileNameFull);
	
	strcpy(MeshNameFull,Path);
	strcat(MeshNameFull,MeshName);
	//	printf("%s\n",MeshNameFull);
	
	strcpy(RstNameFull,Path);
	strcat(RstNameFull,RstName);
	//	printf("%s\n",RstNameFull);
	
	strcpy(OutNameFull,Path);
	strcat(OutNameFull,OutName);
	//	printf("%s\n",OutNameFull);
	
	strcpy(LstNameFull,Path);
	strcat(LstNameFull,LstName);
	//	printf("%s\n",LstNameFull);
	
	strcpy(ConstNameFull,Path);
	strcat(ConstNameFull,ConstName);
	//	printf("%s\n",ConstNameFull);
	
	strcpy(NLNameFull,Path);
	strcat(NLNameFull,NLName);
	//	printf("%s\n",NLNameFull);
	  
	// ----------------------------------------------------------------------
	// set "path" and "Filename" parts of control panel
	// ----------------------------------------------------------------------

	SetDlgItemText(hwnd, IDC_PATH, Path);
	SetDlgItemText(hwnd, IDC_FILENAME, FileName);
	  
	// ----------------------------------------------------------------------
	// Read file.zin and see if it contains restart=YES
	// If so, then file.rst must be present
	// ----------------------------------------------------------------------

	rdu=0;
	fp=fopen(FileNameFull,"r");
	
	if (fp==NULL)
	  MessageBox(hwnd, strcat(Path,FileName), "File not found:", MB_OK | MB_ICONINFORMATION);
	else {

	  while(fgets(string,100,fp) != NULL){
	    i=spliteq(string,left,right);
	    if (i==1){
	      //	      printf("%s %s\n",left,right);
	      if (strcmp(left,"RESTART")==0 && strcmp(right,"YES")==0) rdu=1;
	    }
	  }
	  
	  fclose(fp);
	  Button_Enable(hVfile,TRUE);
	}

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

	if (rdu==0){
	  SetDlgItemText(hwnd, IDC_RSTNAME, "None: initial condition supplied");
	}
	else {
	  fp=fopen(RstNameFull,"r");
	  if (fp==NULL){
	    SetDlgItemText(hwnd, IDC_RSTNAME, "<File not found>");
	    gotfiles=0;
	  }
	  else {
	    fclose(fp);
	    SetDlgItemText(hwnd, IDC_RSTNAME, RstName);
	  }
	}

	fp=fopen(ConstNameFull,"r");
	if (fp==NULL){
	  SetDlgItemText(hwnd, IDC_CONSTNAME, "<No Constants>");
	}
	else {
	  fclose(fp);
	  SetDlgItemText(hwnd, IDC_CONSTNAME, ConstName);
	  Button_Enable(hVconst,TRUE);
	}
	
	fp=fopen(NLNameFull,"r");
	if (fp==NULL){
	  SetDlgItemText(hwnd, IDC_NLNAME, "<No NL file>");
	}
	else {
	  fclose(fp);
	  SetDlgItemText(hwnd, IDC_NLNAME, NLName);
	  Button_Enable(hVNL,TRUE);
	}

	// ----------------------------------------------------------------------
	// set output files text. light up RUN button if we got the files
	// ----------------------------------------------------------------------

	SetDlgItemText(hwnd, IDC_OUTNAME, OutName);
	SetDlgItemText(hwnd, IDC_LSTNAME, LstName);
	
	if (gotfiles==1) 
	  Button_Enable(hRun,TRUE);
	else
	  Button_Enable(hRun,FALSE);
      }
    }
    else if (LOWORD(wParam)==IDC_RUN){
      
      SetDlgItemText(hwnd, IDC_STATUS, "RUNNING");
      Button_Enable(hVlst,FALSE);
      Button_Enable(hRun,FALSE);
      Button_Enable(hOpen,FALSE);
      Button_Enable(hCancel,FALSE);
      
      for (i=0;i<strlen(FileNameFull) && FileNameFull[i]!='.';i++)
	FileNameStem[i]=FileNameFull[i];
      
      FileNameStem[i]='\0';
      
      strcpy(command,"zinc \"");
      strcat(command,FileNameStem);
      strcat(command,"\"");
      
      printf(command);
      
      //      system(command);
      
      params.hwnd=hwnd;
      
      _beginthread(Thread,0,&params);
      
    }
    else if (LOWORD(wParam)==IDC_VFILE)
      LoadTextFileToEdit(hEdit, FileNameFull);

    else if (LOWORD(wParam)==IDC_VMESH)
      LoadTextFileToEdit(hEdit, MeshNameFull);

    else if (LOWORD(wParam)==IDC_VLST)
      LoadTextFileToEdit(hEdit, LstNameFull);

    else if (LOWORD(wParam)==IDC_VCONST)
      LoadTextFileToEdit(hEdit, ConstNameFull);

    else if (LOWORD(wParam)==IDC_VNL)
      LoadTextFileToEdit(hEdit, NLNameFull);

    else if (LOWORD(wParam)==IDC_CANCEL)
      EndDialog(hwnd, 0);

    else
      return FALSE;     // did not handle COMMAND message
  }
  else if (Message==WM_calc_done){
    SetDlgItemText(hwnd, IDC_STATUS, "READY");
    Button_Enable(hVlst,TRUE);
    Button_Enable(hRun,TRUE);
    Button_Enable(hOpen,TRUE);
    Button_Enable(hCancel,TRUE);
  }
  else if (Message==WM_CLOSE)
    EndDialog(hwnd, 0);
  else {
    return FALSE;       // did not handle message
  }

  return TRUE;    // return TRUE if we handled the message
}

// ######################################################################

int spliteq(char* str, char* left, char* right)
{

  // ----------------------------------------------------------------------
  // Split string str at = sign to get null-terminated left and right
  // remove leading blanks from both returned strings. Capitalize both.
  // Return 1 if = found else 0
  // ----------------------------------------------------------------------

  int i,j,ind;

  for (i=0;i<strlen(str);i++){
    if (str[i]=='='){
      ind=0;
      for (j=0;j<i;j++){
	if (str[j]!=' '){
	  left[ind]=toupper(str[j]);
	  ind=ind+1;
	}
      }
      left[ind]='\0';
      
      ind=0;
      for (j=i+1;j<strlen(str);j++){
	if (str[j]!=' ' && str[j]!='\n'){
	  right[ind]=toupper(str[j]);
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

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
	LPSTR lpCmdLine, int nCmdShow)
{
	return DialogBox(hInstance, MAKEINTRESOURCE(IDD_MAIN), NULL, DlgProc);
}

void Thread (PVOID pvoid)
{

  PPARAMS pparams;

  pparams=(PPARAMS) pvoid;

  system(command);

  SendMessage(pparams->hwnd, WM_calc_done, 0, 0);

  _endthread();
}
