#include <iostream>
#include <sys/socket.h>

int main()
{
    //Create socket
    int s = socket(AF_INET, SOCK_STREAM, 0);
    struct sockaddr_in addr {
        AF_INET,
        0x901f, // This has is the reverse of the hex of the port we're gonna use
        0,
    };
    
     
    bind(s, &addr, sizeof(addr));
    //Bind socket to IP port
    //
}
