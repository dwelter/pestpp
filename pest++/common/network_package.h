#ifndef NET_PACKAGE_H_
#define NET_PACKAGE_H_
#include <string>
#include <cstdint>
#include <vector>
#include <memory>

class NetPackage
{
public:
	enum class PackType{UNKN, OK, CONFIRM_OK, READY, REQ_RUNDIR, RUNDIR, REQ_LINPACK, LINPACK, CMD, START_RUN, RUN_FINISHED, RUN_FAILED, RUN_KILLED, TERMINATE,PING,REQ_KILL,IO_ERROR};
	static int get_new_group_id();
	NetPackage(PackType _type=PackType::UNKN, int _group=-1, int _run_id=-1, const std::string &desc="");
	~NetPackage(){}
	const static int DESC_LEN = 41;
	int send(int sockfd, const void *data, unsigned long data_len_l);
	int recv(int sockfd);
	void reset(PackType _type, int _group, int _run_id, const std::string &_desc);
	PackType get_type() const {return type;}
	int get_run_id() const {return run_id;}
	int get_groud_id() const {return group;}
	const std::vector<char> &get_data(){return data;}
	void print_header(std::ostream &fout);
	

private:
	unsigned long data_len;
	static int last_group_id;
	PackType type;
	int group;
	int run_id;
	char desc[DESC_LEN];
	std::vector<char> data;
};

#endif /* NET_PACKAGE_H_ */
