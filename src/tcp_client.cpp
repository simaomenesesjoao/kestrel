#include <sys/socket.h> 
#include <arpa/inet.h> 
#include <iostream>
#include <fstream>
#include <unistd.h> 
#include <string.h> 

#include "tcp_client.hpp"

int tcp_client(std::string *vars){
    // the string 'vars' should contain two special bits in the
    // beginning: '0 1 ...' that indicate the status of the program.
    // The rest is the message to be transmited

    // Status 0 0 - not finished, initializing: not ready to
    //            transmit the information yet
    // Status 1 0 - not finished, initialized: ready to
    //            transmit information about the program
    // Status 1 1 - finished, initialized: program finished 
    //            running. This TCP client should exit soon


    //std::cout << "Entered tcp_client.\n" << std::flush;

    int sock = 0;
    int con, tcp;

    //std::string mes;
    char mes[3000];
    const char *message;

    //char buffer[BUF_SIZE];  //data buffer of 1K  
    //int valread;

    bool connected = false;
    bool finished = false;
    bool str_finished = false;
    bool str_initialized = false;
    while(!finished){

        strcpy(mes, vars->c_str());
        str_finished = mes[2]=='1';
        str_initialized = mes[0]=='1';

        // Attempt to connect continuously
        if(!connected ){
            if(!str_finished){
                con = connect(&sock);
                if(con == 0){
                    connected = true;
                    std::cout << "Connection established.\n";
                }
            } else {

                finished = true;
            }
        }

        if(connected){
            if(str_finished){
                tcp = send_tcp(sock, (char*)mes);

                close(sock);
                finished = true;
                continue;
            }
            if(str_initialized){
                tcp = send_tcp(sock, (char*)mes);
                if(tcp == -1){
                    connected = false;
                    std::cout << "Connection lost.\n";
                    close(sock);
                }
            }

            if(!str_initialized){
                //char *end_message = "Initializing.";
                //tcp = send_tcp(sock, end_message);
                tcp = send_tcp(sock, (char*)mes);
                //valread = read( sock , buffer, 1024);

                //close(sock);
                //finished = true;
                //continue;
            }

        }

        for(int s = 0; s < 5; s++){
            strcpy(mes, vars->c_str());
            str_finished = mes[2]=='1';
            sleep(100);
            if(!connected or str_finished)
                break;
        }
    }

    //std::cout << "Left tcp_client.\n" << std::flush;
    return 0; 
} 

int connect(int *sock){

    unsigned PORT = 60501;
    std::string IP ="127.0.0.1";

    char *home;
    home = getenv("HOME");
    if(home != NULL){
        std::string PORT_str;
        std::ifstream file;
        file.open(std::string(home) + "/.config/kestrel/ip");
        if(!file.fail()){
            std::getline(file, IP);
        }
        file.close();
        file.open(std::string(home) + "/.config/kestrel/port");
        if(!file.fail()){
            std::getline(file, PORT_str);
            PORT = std::stoi(PORT_str);
        }
        file.close();
    }


    struct sockaddr_in serv_addr; 
    if ((*sock = socket(AF_INET, SOCK_STREAM, 0)) < 0){ 
        printf("\n Socket creation error \n"); 
        return -1; 
    } 
   
    serv_addr.sin_family = AF_INET; 
    serv_addr.sin_port = htons(PORT); 
       
    // Convert IPv4 and IPv6 addresses from text to binary form 
    if(inet_pton(AF_INET, (char*)IP.c_str(), &serv_addr.sin_addr)<=0){ 
        printf("\nInvalid address/ Address not supported \n"); 
        return -1; 
    } 
   

    if (connect(*sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0){ 
        std::cout << "Connection to " << IP << ":" << PORT << " failed.\n";
        return -1; 
    } 
    return 0;
}

int send_tcp(int socket, char *message){
    //std::cout << "Entered send_tcp.\n";
    int send_result = send(socket, message, strlen(message), MSG_NOSIGNAL); 
    //std::cout << "Left send_tcp.\n";
    return send_result;
}
