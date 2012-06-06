/*
 Copyright 2012, S.S. Papadopulos & Assoc. Inc.

    This file is part of GENIE.

    The GENIE is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GENIE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GENIE.  If not, see<http://www.gnu.org/licenses/>.
*/

// cmuffels, Genie version 1, 2012

#include <iostream>
#include <time.h>

#include "message.h"

using std::cout;

//****************************************************************************
void MESSAGE::alloc()
//****************************************************************************
//CTM Jan 2011

/*
      Allocate memory space for data portion of message
*/

{
  data=new BUFFER;
  data->alloc(datasize);
}

//****************************************************************************
void MESSAGE::dealloc()
//****************************************************************************
//CTM Jan 2011

/*
      Routine to clear data memory block
*/

{
  //cout<<" deallocating message data..."<<'\n';
  data->dealloc();
  //cout<<" deleting message data..."<<'\n';
  delete data;
  //cout<<" clearing data pointer..."<<'\n';
  data=NULL;

}

//****************************************************************************
void MESSAGE::setmsgsize(size_t _size)
//****************************************************************************
//CTM Jan 2011

/*
      Routine to set message size
*/

{

  msgsize=_size;

}

//****************************************************************************
void MESSAGE::setdatasize(size_t _size)
//****************************************************************************
//CTM Jan 2011

/*
      Routine to set size of data portion of message
*/

{

  datasize=_size;

}

//****************************************************************************
int MESSAGE::receiveheader(SOCKET &_socket)
//****************************************************************************
//CTM Jan 2011

/*
      Receive bytes in header
*/

{
  int nrecv,nleft,ntry;
  bool haveheader=false;
  BUFFER buffer;
  fd_set fdset;

  // assign socket to fdset
  FD_ZERO(&fdset);
  FD_SET(_socket,&fdset);

  // initialize header buffer
  buffer.alloc(sizeof(header));

  try
  {
    // initialize byte counters
    nrecv=0;
    nleft=sizeof(header);
    while(1)
    {

      if(haveheader) break;

      // wait for message to come in on socket
      switch(waitformessage(_socket,COMM_WAIT_TIME,0))
      {
      case SELECT_TIMEOUT:
        return MSG_NONE;
        break;
      case SOCKET_ERROR:
        throw -1;
        break;
      default:
        //receive some bytes
        ntry=recv(_socket,&buffer.buf[nrecv],nleft,0);
        //cout<<"\nNTRY: "<<ntry<<'\n';
        switch(ntry)
        {
        case SOCKET_ERROR:
          throw -2;
          break;
        case CONNECTION_LOST:
          throw -3;
          break;
        default:
          nrecv+=ntry;
          nleft-=ntry;
          if(nleft==0)
          {
            //copy buf to header
            buffer.copyto(0,(char*)&header,sizeof(header));
            buffer.dealloc();
            setdatasize(header.bytesize*header.nbytes);
            setmsgsize(sizeof(header)+datasize);
            //all done
            haveheader=true;
          }
          break;
        } //end of bytes received switch
        break;
      }  // end of waitformessage switch
    } // end of while loop

    return 1;
  }
  catch (int err)
  {
    if(err==-1)
      cout<<"\n\n Unexpected SOCKET error waiting for message header."<<'\n';
    else if(err==-2)
      cout<<"\n\n Unexpected SOCKET error receiving message header."<<'\n';
    else if(err==-3)
      cout<<"\n\n Client disconnected."<<'\n';
    else if(err==-10)
      cout<<"\n\n Unexpected SOCKET time out waiting for message header."<<'\n';
    buffer.dealloc();
    return 0;
  }
}

//****************************************************************************
int MESSAGE::receivedata(SOCKET &_socket)
//****************************************************************************
//CTM Jan 2011

/*
      Routine description here
*/

{

  size_t nrecv,nleft,ntry;
  bool havedata=false;

  try
  {
    //initialize byte counters
    nrecv=0;
    nleft=datasize;
    while(1)
    {

      if(havedata) break;

      // wait for message to come in on socket
      switch(waitformessage(_socket,-1,0))
      {
      case SELECT_TIMEOUT:
        throw -10;
        break;
      case SOCKET_ERROR:
        throw -1;
        break;
      default:
        //receive some bytes
        ntry=(size_t)recv(_socket,&data->buf[nrecv],(int)nleft,0);
        //cout<<nrecv<<","<<nleft<<","<<ntry<<'\n';
        switch(ntry)
        {
        case SOCKET_ERROR:
          throw -2;
          break;
        case CONNECTION_LOST:
          throw -3;
          break;
        default:
          nrecv+=ntry;
          nleft-=ntry;
          if(nleft==0) havedata=true;
          break; //all done.
        }  // end of bytes received switch
        break;
      }  // end of waitformessage switch
    }  // end of while
    return 1;
  }  // end of try
  catch (int err)
  {
    if(err==-1)
      cout<<"\n\n Unexpected SOCKET error waiting for message data."<<'\n';
    if(err==-2)
      cout<<"\n\n Unexpected SOCKET error receiving message data."<<'\n';
    if(err==-3)
    {
      cout<<"\n\n Client disconnected unexpectedly while";
      cout<<" receiving message data."<<'\n';
    }
    if(err==-10)
      cout<<"\n\n Unexpected SOCKET time out waiting for message data."<<'\n';
    return 0;
  }
}

//****************************************************************************
int MESSAGE::sendme(SOCKET &_socket)
//****************************************************************************
//CTM Jan 2011

/*
      Send message to socket
*/

{

  size_t nleft,nsent,ntry;
  BUFFER buffer;

  // initialize buffer to contain both header and data
  buffer.alloc(msgsize);
  buffer.copyfrom(0,(char*)&header,sizeof(header));
  buffer.copyfrom(sizeof(header),data->buf,datasize);

  try
  {
    nleft=msgsize;
    nsent=0;
    while(nleft>0)
    {
      ntry=(size_t)send(_socket,&buffer.buf[nsent],(int)nleft,0);
      if(ntry==SOCKET_ERROR) throw -1;
      nleft-=ntry;
      nsent+=ntry;
      //cout<<nsent<<","<<nleft<<","<<ntry<<'\n';
    }
    buffer.dealloc();
    return 1;
  }
  catch(...)
  {
    cout<<"\n\n Error sending buffer. Unexpected SOCKET error. "<<_socket<<'\n';
    buffer.dealloc();
    return 0;
  }

}

//****************************************************************************
int MESSAGE::waitformessage(SOCKET &_socket,long sec,long usec)
//****************************************************************************
//CTM Dec 2010
// from http://www.winsocketdotnetworkprogramming.com/winsock2programming/
//                                             winsock2advancedcode1c.html

/*
      Routine to wait sec and usec long for activity on a socket
*/
{

  // Setup timeval variable
  struct timeval timeout;
  fd_set fdset;

  // assign socket to fdset
  FD_ZERO(&fdset);
  FD_SET(_socket,&fdset);

  // assign the second and microsecond variables
  timeout.tv_sec = sec;
  timeout.tv_usec = usec;

  // Possible return values:
  // -1: error occurred
  // 0: timed out
  // > 0: data ready to be read
  if(sec==-1)
    return select(FD_SETSIZE, &fdset, NULL, NULL, NULL);
  else
    return select(FD_SETSIZE, &fdset, NULL, NULL, &timeout);

}